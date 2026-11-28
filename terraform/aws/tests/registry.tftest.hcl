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
    condition = jsondecode(aws_iam_role_policy.registry_auth[0].policy).Statement == [{
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
    condition     = jsondecode(aws_iam_role_policy.registry_auth[0].policy).Statement[1].Resource == [var.ghcr_secret_kms_key_arn] && jsondecode(aws_iam_role_policy.registry_auth[0].policy).Statement[1].Condition.StringEquals["kms:EncryptionContext:SecretARN"] == var.ghcr_secret_arn
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

run "public_images_without_registry_credentials" {
  command = plan
  variables {
    ghcr_secret_arn         = null
    basecall_instance_types = ["g6e"]
    gpu_instance_types      = ["g6e"]
    bam_queue_enabled       = true
  }
  assert {
    condition     = length(aws_iam_role_policy.registry_auth) == 0
    error_message = "Anonymous pulls must not grant access to registry secrets."
  }
  assert {
    condition = alltrue([for user_data in [aws_launch_template.batch.user_data, aws_launch_template.bam[0].user_data, aws_launch_template.bam_low[0].user_data] :
      !strcontains(base64decode(user_data), "20-ghcr-auth.conf") && !strcontains(base64decode(user_data), "--secret-arn") && strcontains(base64decode(user_data), "systemctl start ecs")
    ])
    error_message = "All worker templates must start ECS without registry authentication for public images."
  }
}

run "reject_registry_key_without_secret" {
  command = plan
  variables {
    ghcr_secret_arn         = null
    ghcr_secret_kms_key_arn = "arn:aws:kms:us-east-1:123456789012:key/test-key"
  }
  expect_failures = [var.ghcr_secret_kms_key_arn]
}

run "bundled_cli_runtime_library" {
  command = plan
  variables {
    ghcr_secret_arn = null
    bam_queue_enabled = true
  }
  assert {
    condition = alltrue([for user_data in [aws_launch_template.batch.user_data, aws_launch_template.bam[0].user_data, aws_launch_template.bam_low[0].user_data] :
      strcontains(base64decode(user_data), "dnf install -y unzip zlib") &&
      strcontains(base64decode(user_data), "cp -L /lib64/libz.so.1 /opt/aws-cli/v2/current/runtime-libs/libz.so.1") &&
      strcontains(base64decode(user_data), "export LD_LIBRARY_PATH=") &&
      strcontains(base64decode(user_data), "exec \"$cli_root/dist/aws\" \"$@\"")
    ])
    error_message = "All downloaded CLI installations must bundle zlib and invoke the original executable through the library-path wrapper."
  }
}
