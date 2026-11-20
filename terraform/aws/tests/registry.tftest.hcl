mock_provider "aws" {
  mock_data "aws_s3_bucket" {
    defaults = { id = "test-bucket", arn = "arn:aws:s3:::test-bucket" }
  }
  mock_data "aws_partition" {
    defaults = { partition = "aws", dns_suffix = "amazonaws.com" }
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
  bucket_name     = "test-bucket"
  ghcr_secret_arn = "arn:aws:secretsmanager:us-east-1:123456789012:secret:test-ghcr-Ab1234"
}
run "registry_worker_policy_and_bootstrap" {
  command = plan
  variables { gpu_instance_types = ["g6e"] }
  assert {
    condition = jsondecode(aws_iam_role_policy.registry_auth.policy).Statement == [{
      Effect = "Allow", Action = ["secretsmanager:GetSecretValue"], Resource = [var.ghcr_secret_arn]
    }]
    error_message = "Default registry access must be limited to one existing secret."
  }
  assert {
    condition     = strcontains(base64decode(aws_launch_template.batch.user_data), "ExecStartPre=/usr/bin/python3 /opt/oncoseq-leuk-registry/configure.py --aws-cli /opt/aws-cli/v2/current/bin/aws --region us-east-1 --secret-arn ${var.ghcr_secret_arn}")
    error_message = "Workers must fetch the credential before ECS starts."
  }
}
run "registry_customer_key_and_existing_cli" {
  command = plan
  variables {
    ghcr_secret_kms_key_arn = "arn:aws:kms:us-east-1:123456789012:key/test-key"
    use_existing_aws_cli    = true
  }
  assert {
    condition     = jsondecode(aws_iam_role_policy.registry_auth.policy).Statement[1].Resource == [var.ghcr_secret_kms_key_arn] && jsondecode(aws_iam_role_policy.registry_auth.policy).Statement[1].Condition.StringEquals["kms:EncryptionContext:SecretARN"] == var.ghcr_secret_arn
    error_message = "Decrypt must be scoped to the existing secret key and encryption context."
  }
  assert {
    condition     = strcontains(base64decode(aws_launch_template.batch.user_data), "--aws-cli ${var.existing_aws_cli_path}")
    error_message = "Registry auth must honor the existing CLI installation."
  }
}
run "reject_empty_secret" {
  command = plan
  variables { ghcr_secret_arn = "" }
  expect_failures = [var.ghcr_secret_arn]
}
run "reject_wrong_secret_region" {
  command = plan
  variables { ghcr_secret_arn = "arn:aws:secretsmanager:us-west-2:123456789012:secret:test-ghcr-Ab1234" }
  expect_failures = [var.ghcr_secret_arn]
}
