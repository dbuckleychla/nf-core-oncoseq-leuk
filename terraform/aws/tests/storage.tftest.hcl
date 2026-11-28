# Offline plans: mocked AWS provider; no AWS credentials or deployments.
mock_provider "aws" {
  mock_data "aws_s3_bucket" {
    defaults = {
      id  = "existing-test-bucket"
      arn = "arn:aws:s3:::existing-test-bucket"
    }
  }
  mock_data "aws_partition" {
    defaults = { partition = "aws" }
  }
  mock_data "aws_caller_identity" {
    defaults = { account_id = "123456789012" }
  }
  mock_data "aws_availability_zones" {
    defaults = { names = ["us-east-1a", "us-east-1b"] }
  }
  mock_data "aws_iam_policy_document" {
    defaults = { json = "{}" }
  }
}
variables {
  cpu_ami_id      = "ami-0123456789abcdef0"
  gpu_ami_id      = "ami-0123456789abcdef1"
  region          = "us-east-1"
  bucket_name     = "existing-test-bucket"
  ghcr_secret_arn = "arn:aws:secretsmanager:us-east-1:123456789012:secret:test-ghcr-Ab1234"
}
run "default_prefix" {
  command = plan
  # Deployment tfvars are loaded by terraform test; isolate this fixture's prefix.
  variables { s3_prefix = "oncoseq-leuk/" }
  assert {
    condition     = output.work_dir == "s3://existing-test-bucket/oncoseq-leuk/work"
    error_message = "Work must use the required existing bucket and case-sensitive default prefix."
  }
  assert {
    condition     = output.nextflow_params.outdir == "s3://existing-test-bucket/oncoseq-leuk/results"
    error_message = "Results must use the same default prefix."
  }
}
run "custom_prefix" {
  command = plan
  variables { s3_prefix = "team/run-1/" }
  assert {
    condition     = output.work_dir == "s3://existing-test-bucket/team/run-1/work"
    error_message = "Normalize trailing slash without duplicating separators."
  }
  assert {
    condition     = output.nextflow_params.outdir == "s3://existing-test-bucket/team/run-1/results"
    error_message = "Custom prefix must apply to results too."
  }
}
run "reject_root_prefix" {
  command = plan
  variables { s3_prefix = "/" }
  expect_failures = [var.s3_prefix]
}
run "reject_parent_prefix" {
  command = plan
  variables { s3_prefix = "team/../other" }
  expect_failures = [var.s3_prefix]
}
run "reject_bucket_url" {
  command = plan
  variables { bucket_name = "s3://existing-test-bucket" }
  expect_failures = [var.bucket_name]
}
run "worker_diagnostics_and_storage" {
  command = plan
  variables { gpu_instance_types = ["g6e"] }
  assert {
    condition     = aws_launch_template.batch.block_device_mappings[0].ebs[0].iops == 6000 && aws_launch_template.batch.block_device_mappings[0].ebs[0].throughput == 500
    error_message = "CPU/GPU shared launch template must provision explicit gp3 performance."
  }
  assert {
    condition     = aws_cloudwatch_log_group.hosts.retention_in_days == var.log_retention_days
    error_message = "Host logs must retain diagnostics after worker termination."
  }
}
run "reject_disk_ratio" {
  command = plan
  variables {
    scratch_iops             = 3000
    scratch_throughput_mibps = 1000
  }
  expect_failures = [var.scratch_throughput_mibps]
}

run "exclusive_bam_workers" {
  command = plan
  override_resource {
    target          = aws_batch_compute_environment.bam[0]
    override_during = plan
    values          = { arn = "arn:aws:batch:us-east-1:123456789012:compute-environment/bam-high" }
  }
  override_resource {
    target          = aws_batch_compute_environment.bam_low[0]
    override_during = plan
    values          = { arn = "arn:aws:batch:us-east-1:123456789012:compute-environment/bam-low" }
  }
  variables {
    bam_queue_enabled  = true
    bam_max_vcpus      = 256
    bam_scratch_gb     = 2000
    bam_low_max_vcpus  = 64
    bam_low_scratch_gb = 2000
    scratch_gb         = 500
  }
  assert {
    condition     = length(aws_batch_compute_environment.bam[0].compute_resources[0].instance_type) == 1 && contains(aws_batch_compute_environment.bam[0].compute_resources[0].instance_type, "r8i.8xlarge")
    error_message = "BAM queue must use only the fixed 32-vCPU instance type."
  }
  assert {
    condition     = aws_launch_template.bam[0].block_device_mappings[0].ebs[0].volume_size == 2000
    error_message = "BAM workers must have separate 2000 GiB staging disks."
  }
  assert {
    condition     = output.nextflow_params.aws_bam_high_vcpus == "32" && output.nextflow_params.aws_bam_high_queue == "${var.name}-bam"
    error_message = "Retain the existing high queue name and export its full-instance reservation."
  }
  assert {
    condition     = length(aws_batch_compute_environment.bam_low[0].compute_resources[0].instance_type) == 1 && contains(aws_batch_compute_environment.bam_low[0].compute_resources[0].instance_type, "r8i.2xlarge")
    error_message = "BAM low must use only fixed 8-vCPU r8i.2xlarge workers."
  }
  assert {
    condition     = aws_launch_template.bam_low[0].block_device_mappings[0].ebs[0].volume_size == 2000 && aws_launch_template.bam_low[0].block_device_mappings[0].ebs[0].iops == var.scratch_iops && aws_launch_template.bam_low[0].block_device_mappings[0].ebs[0].throughput == var.scratch_throughput_mibps
    error_message = "Low workers must keep 2000-GiB staging disks and shared gp3 performance."
  }
  assert {
    condition     = aws_launch_template.bam_low[0].user_data == aws_launch_template.bam[0].user_data && aws_launch_template.bam_low[0].metadata_options[0].http_tokens == "required"
    error_message = "Both fleets must retain registry/bootstrap diagnostics and IMDSv2."
  }
  assert {
    condition     = output.nextflow_params.aws_bam_low_vcpus == "8" && output.nextflow_params.aws_bam_low_queue == "${var.name}-bam-low" && output.nextflow_params.aws_bam_low_queue != output.nextflow_params.aws_bam_high_queue
    error_message = "Export distinct high/low queues with full 32/8-vCPU reservations."
  }
  assert {
    condition     = length(aws_batch_job_queue.bam[0].compute_environment_order) == 1 && one(aws_batch_job_queue.bam[0].compute_environment_order).compute_environment == aws_batch_compute_environment.bam[0].arn && length(aws_batch_job_queue.bam_low[0].compute_environment_order) == 1 && one(aws_batch_job_queue.bam_low[0].compute_environment_order).compute_environment == aws_batch_compute_environment.bam_low[0].arn
    error_message = "Each BAM queue must use only its own fleet to preserve exclusive placement."
  }
  assert {
    condition     = aws_batch_compute_environment.bam_low[0].compute_resources[0].min_vcpus == 0 && aws_batch_compute_environment.bam_low[0].compute_resources[0].max_vcpus == var.bam_low_max_vcpus && aws_batch_compute_environment.bam_low[0].compute_resources[0].type == "EC2"
    error_message = "Low capacity must be independently configured, on demand and scale to zero."
  }
  assert {
    condition     = aws_launch_template.batch.block_device_mappings[0].ebs[0].volume_size == 500
    error_message = "Dedicated disk changes must not resize existing CPU/GPU workers."
  }
}

run "bam_queues_disabled" {
  command = plan
  variables { bam_queue_enabled = false }
  assert {
    condition     = length(aws_batch_compute_environment.bam) == 0 && length(aws_batch_compute_environment.bam_low) == 0 && output.nextflow_params.aws_bam_high_queue == null && output.nextflow_params.aws_bam_low_queue == null
    error_message = "Disabling the optional BAM fleets must retain the original general/GPU routing."
  }
}
run "custom_bam_low_storage" {
  command = plan
  variables {
    bam_queue_enabled  = true
    bam_low_max_vcpus  = 16
    bam_low_scratch_gb = 3000
    bam_scratch_gb     = 2500
  }
  assert {
    condition     = aws_launch_template.bam_low[0].block_device_mappings[0].ebs[0].volume_size == 3000 && aws_launch_template.bam[0].block_device_mappings[0].ebs[0].volume_size == 2500 && aws_batch_compute_environment.bam_low[0].compute_resources[0].max_vcpus == 16
    error_message = "Low storage/capacity must be independent of existing high settings."
  }
}
run "reject_low_capacity" {
  command = plan
  variables { bam_low_max_vcpus = 12 }
  expect_failures = [var.bam_low_max_vcpus]
}
run "reject_low_disk" {
  command = plan
  variables { bam_low_scratch_gb = 499 }
  expect_failures = [var.bam_low_scratch_gb]
}

run "independent_gpu_queues" {
  command = plan
  variables {
    basecall_instance_types = ["g6e"]
    gpu_instance_types      = ["g5"]
    basecall_max_vcpus      = 64
    gpu_max_vcpus           = 32
  }
  assert {
    condition     = output.nextflow_params.aws_basecall_queue == "${var.name}-basecall" && output.nextflow_params.aws_gpu_queue == "${var.name}-gpu"
    error_message = "Basecalling and general GPU queues must have distinct exported names."
  }
  assert {
    condition     = aws_batch_compute_environment.basecall[0].compute_resources[0].max_vcpus == 64 && aws_batch_compute_environment.gpu[0].compute_resources[0].max_vcpus == 32
    error_message = "GPU fleets must have independently configurable capacities."
  }
  assert {
    condition     = aws_batch_compute_environment.basecall[0].compute_resources[0].ec2_configuration[0].image_id_override == var.gpu_ami_id && aws_batch_compute_environment.this.compute_resources[0].ec2_configuration[0].image_id_override == var.cpu_ami_id
    error_message = "Compute environments must use the supplied pinned AMIs."
  }
}

run "basecalling_without_general_gpu_workers" {
  command = plan
  variables {
    basecall_instance_types = ["g6e"]
    gpu_instance_types      = []
  }
  assert {
    condition     = output.nextflow_params.aws_basecall_queue != null && output.nextflow_params.aws_gpu_queue == null
    error_message = "Basecalling capacity must be independently selectable."
  }
}
