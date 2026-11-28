# BAM high: retain existing resource addresses/names for state-compatible rollout.
# Each task reserves all 32 vCPUs; the low fleet below reserves all eight.
resource "aws_launch_template" "bam" {
  count       = var.bam_queue_enabled ? 1 : 0
  name_prefix = "${var.name}-bam-"
  user_data = base64encode(templatefile("${path.module}/templates/user_data.mime", {
    aws_cli_version       = var.aws_cli_version
    use_existing_aws_cli  = var.use_existing_aws_cli
    existing_aws_cli_path = var.existing_aws_cli_path
    host_agent_config     = local.host_agent_config
    region                = var.region
    ghcr_secret_arn       = var.ghcr_secret_arn
    registry_auth_script  = file("${path.module}/templates/configure_registry_auth.py")
  }))
  block_device_mappings {
    device_name = "/dev/xvda"
    ebs {
      volume_size           = var.bam_scratch_gb
      volume_type           = "gp3"
      iops                  = var.scratch_iops
      throughput            = var.scratch_throughput_mibps
      encrypted             = true
      delete_on_termination = true
    }
  }
  metadata_options {
    http_tokens                 = "required"
    http_put_response_hop_limit = 2
  }
}
resource "aws_batch_compute_environment" "bam" {
  count        = var.bam_queue_enabled ? 1 : 0
  name_prefix  = "${var.name}-bam-"
  type         = "MANAGED"
  service_role = aws_iam_role.batch.arn
  compute_resources {
    type                = "EC2"
    allocation_strategy = "BEST_FIT_PROGRESSIVE"
    min_vcpus           = 0
    max_vcpus           = var.bam_max_vcpus
    instance_type       = ["r8i.8xlarge"]
    instance_role       = aws_iam_instance_profile.batch.arn
    spot_iam_fleet_role = null
    subnets             = local.subnets
    security_group_ids  = [aws_security_group.batch.id]
    ec2_configuration {
      image_type        = "ECS_AL2023"
      image_id_override = var.cpu_ami_id
    }
    launch_template {
      launch_template_id = aws_launch_template.bam[0].id
      version            = tostring(aws_launch_template.bam[0].latest_version)
    }
  }
  depends_on = [aws_iam_role_policy.registry_auth, aws_iam_role_policy.host_diagnostics, aws_iam_role_policy_attachment.batch, aws_iam_role_policy_attachment.instance, aws_iam_role_policy_attachment.spot]
  lifecycle { create_before_destroy = true }
}
resource "aws_batch_job_queue" "bam" {
  count    = var.bam_queue_enabled ? 1 : 0
  name     = "${var.name}-bam"
  state    = "ENABLED"
  priority = 1
  compute_environment_order {
    order               = 1
    compute_environment = aws_batch_compute_environment.bam[0].arn
  }
}

# BAM low: a separate fixed-size fleet with no shared compute environments.
resource "aws_launch_template" "bam_low" {
  count       = var.bam_queue_enabled ? 1 : 0
  name_prefix = "${var.name}-bam-low-"
  user_data = base64encode(templatefile("${path.module}/templates/user_data.mime", {
    aws_cli_version       = var.aws_cli_version
    use_existing_aws_cli  = var.use_existing_aws_cli
    existing_aws_cli_path = var.existing_aws_cli_path
    host_agent_config     = local.host_agent_config
    region                = var.region
    ghcr_secret_arn       = var.ghcr_secret_arn
    registry_auth_script  = file("${path.module}/templates/configure_registry_auth.py")
  }))
  block_device_mappings {
    device_name = "/dev/xvda"
    ebs {
      volume_size           = var.bam_low_scratch_gb
      volume_type           = "gp3"
      iops                  = var.scratch_iops
      throughput            = var.scratch_throughput_mibps
      encrypted             = true
      delete_on_termination = true
    }
  }
  metadata_options {
    http_tokens                 = "required"
    http_put_response_hop_limit = 2
  }
}
resource "aws_batch_compute_environment" "bam_low" {
  count        = var.bam_queue_enabled ? 1 : 0
  name_prefix  = "${var.name}-bam-low-"
  type         = "MANAGED"
  service_role = aws_iam_role.batch.arn
  compute_resources {
    type                = "EC2"
    allocation_strategy = "BEST_FIT_PROGRESSIVE"
    min_vcpus           = 0
    max_vcpus           = var.bam_low_max_vcpus
    instance_type       = ["r8i.2xlarge"]
    instance_role       = aws_iam_instance_profile.batch.arn
    spot_iam_fleet_role = null
    subnets             = local.subnets
    security_group_ids  = [aws_security_group.batch.id]
    ec2_configuration {
      image_type        = "ECS_AL2023"
      image_id_override = var.cpu_ami_id
    }
    launch_template {
      launch_template_id = aws_launch_template.bam_low[0].id
      version            = tostring(aws_launch_template.bam_low[0].latest_version)
    }
  }
  depends_on = [aws_iam_role_policy.registry_auth, aws_iam_role_policy.host_diagnostics, aws_iam_role_policy_attachment.batch, aws_iam_role_policy_attachment.instance, aws_iam_role_policy_attachment.spot]
  lifecycle { create_before_destroy = true }
}
resource "aws_batch_job_queue" "bam_low" {
  count    = var.bam_queue_enabled ? 1 : 0
  name     = "${var.name}-bam-low"
  state    = "ENABLED"
  priority = 1
  compute_environment_order {
    order               = 1
    compute_environment = aws_batch_compute_environment.bam_low[0].arn
  }
}
