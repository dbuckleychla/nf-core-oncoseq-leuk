output "nextflow_params" {
  value = {
    aws_bam_high_queue = var.bam_queue_enabled ? aws_batch_job_queue.bam[0].name : null
    aws_bam_high_vcpus = var.bam_queue_enabled ? "32" : null
    aws_bam_low_queue  = var.bam_queue_enabled ? aws_batch_job_queue.bam_low[0].name : null
    aws_bam_low_vcpus  = var.bam_queue_enabled ? "8" : null
    aws_basecall_queue = length(aws_batch_job_queue.basecall) > 0 ? aws_batch_job_queue.basecall[0].name : null
    aws_gpu_queue      = length(aws_batch_job_queue.gpu) > 0 ? aws_batch_job_queue.gpu[0].name : null
    aws_region         = var.region
    aws_queue          = aws_batch_job_queue.this.name
    aws_job_role       = aws_iam_role.task.arn
    aws_logs_group     = aws_cloudwatch_log_group.batch.name
    aws_cli_path       = var.use_existing_aws_cli ? var.existing_aws_cli_path : "/opt/aws-cli/v2/current/bin/aws"
    outdir             = "s3://${local.bucket_name}/${local.results_prefix}"
  }
}
output "work_dir" { value = "s3://${local.bucket_name}/${local.work_prefix}" }
output "coordinator_policy_arn" { value = aws_iam_policy.coordinator.arn }
output "vpc_id" { value = local.vpc_id }
output "subnet_ids" { value = local.subnets }
