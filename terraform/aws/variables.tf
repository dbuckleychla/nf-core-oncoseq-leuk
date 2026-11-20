variable "region" {
  type        = string
  description = "AWS region. Credentials use the standard AWS SDK chain."
}
variable "ghcr_secret_arn" {
  type        = string
  nullable    = false
  description = "Required existing Secrets Manager ARN for GHCR: JSON username plus accessToken (MSQ format) or password. Supply in tfvars or TF_VAR_ghcr_secret_arn; Terraform prompts if omitted. Never supply the token itself."
  validation {
    condition     = can(regex("^arn:aws[a-z-]*:secretsmanager:[a-z0-9-]+:[0-9]{12}:secret:[A-Za-z0-9/_+=.@-]+$", var.ghcr_secret_arn)) && try(split(":", var.ghcr_secret_arn)[3] == var.region, false)
    error_message = "Supply a Secrets Manager secret ARN in the deployment region."
  }
}
variable "ghcr_secret_kms_key_arn" {
  type        = string
  default     = null
  description = "Existing customer-managed KMS key encrypting the GHCR secret, if used. Null for the default Secrets Manager key. Grants worker decrypt only."
  validation {
    condition     = var.ghcr_secret_kms_key_arn == null ? true : can(regex("^arn:aws[a-z-]*:kms:[a-z0-9-]+:[0-9]{12}:key/[A-Za-z0-9-]+$", var.ghcr_secret_kms_key_arn))
    error_message = "Use a KMS key ARN, not an alias, or null."
  }
}
variable "name" {
  type    = string
  default = "oncoseq-leuk"
  validation {
    condition     = can(regex("^[a-z][a-z0-9-]{2,23}$", var.name))
    error_message = "Use 3–24 lower-case letters, digits and hyphens."
  }
}
variable "tags" {
  type    = map(string)
  default = {}
}
variable "vpc_id" {
  type        = string
  default     = null
  description = "Existing VPC; null creates a VPC with private Batch subnets and NAT."
}
variable "subnet_ids" {
  type    = list(string)
  default = []
  validation {
    condition     = (var.vpc_id == null && length(var.subnet_ids) == 0) || (var.vpc_id != null && length(var.subnet_ids) > 0)
    error_message = "Supply both vpc_id and subnet_ids, or neither."
  }
}
variable "vpc_cidr" {
  type    = string
  default = "10.83.0.0/16"
}
variable "bucket_name" {
  type        = string
  nullable    = false
  description = "Required existing S3 bucket name. Terraform never creates, deletes or configures this bucket."
  validation {
    condition     = can(regex("^[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]$", var.bucket_name))
    error_message = "Supply an existing S3 bucket name, not an ARN, URL or path."
  }
}
variable "s3_prefix" {
  type        = string
  default     = "oncoseq-leuk/"
  nullable    = false
  description = "Nonempty relative prefix containing work/ and results/; optional trailing slash."
  validation {
    condition     = can(regex("^[A-Za-z0-9][A-Za-z0-9._/-]*/?$", var.s3_prefix)) && alltrue([for part in split("/", trimsuffix(var.s3_prefix, "/")) : !contains(["", ".", ".."], part)])
    error_message = "Use a nonempty relative S3 prefix without empty, dot or parent-directory components."
  }
}
variable "read_bucket_arns" {
  type        = list(string)
  default     = []
  description = "Additional input/reference bucket ARNs readable by tasks and coordinator."
}
variable "kms_key_arn" {
  type        = string
  default     = null
  description = "Optional existing key used by the supplied S3 bucket. Grants task/coordinator KMS access only; does not configure encryption."
}
variable "read_kms_key_arns" {
  type    = list(string)
  default = []
}
variable "permissions_boundary_arn" {
  type    = string
  default = null
}
variable "instance_types" {
  type    = list(string)
  default = ["m6i", "r6i"]
}
variable "max_vcpus" {
  type    = number
  default = 128
  validation {
    condition     = var.max_vcpus > 0
    error_message = "max_vcpus must be positive."
  }
}
variable "spot" {
  type    = bool
  default = false
}
variable "scratch_gb" {
  type    = number
  default = 500
}
variable "log_retention_days" {
  type    = number
  default = 30
}
variable "aws_cli_version" {
  type    = string
  default = "2.31.0"
  validation {
    condition     = can(regex("^2\\.[0-9]+\\.[0-9]+$", var.aws_cli_version))
    error_message = "Use an explicit AWS CLI v2 release version."
  }
}
variable "gpu_instance_types" {
  type        = list(string)
  default     = []
  description = "Optional x86-64 NVIDIA instance families, e.g. [g5]; empty disables GPU resources."
}
variable "gpu_max_vcpus" {
  type    = number
  default = 32
}

variable "use_existing_aws_cli" {
  type        = bool
  default     = false
  description = "Use AWS CLI already installed in the EC2 AMI instead of downloading it during bootstrap. Applies to CPU and GPU instances."
}
variable "existing_aws_cli_path" {
  type        = string
  default     = "/usr/local/aws-cli/v2/current/bin/aws"
  description = "Absolute path to the existing, self-contained AWS CLI installation in both EC2 AMIs; Nextflow mounts this installation into task containers."
  validation {
    condition     = can(regex("^/[A-Za-z0-9_./-]+/aws$", var.existing_aws_cli_path))
    error_message = "Use an absolute AWS CLI executable path ending in /aws, without spaces or shell metacharacters."
  }
}

variable "scratch_iops" {
  type        = number
  default     = 6000
  description = "Provisioned gp3 root-disk IOPS shared by Docker and task scratch; incurs EBS charges above baseline."
  validation {
    condition     = var.scratch_iops >= 3000 && var.scratch_iops <= 16000 && floor(var.scratch_iops) == var.scratch_iops
    error_message = "scratch_iops must be an integer from 3000 to 16000."
  }
}
variable "scratch_throughput_mibps" {
  type        = number
  default     = 500
  description = "gp3 throughput in MiB/s shared by image extraction and task I/O; limited by the EC2 instance EBS bandwidth."
  validation {
    condition     = var.scratch_throughput_mibps >= 125 && var.scratch_throughput_mibps <= 1000 && var.scratch_throughput_mibps <= var.scratch_iops / 4 && floor(var.scratch_throughput_mibps) == var.scratch_throughput_mibps
    error_message = "Use integer throughput 125–1000 MiB/s, at most one quarter of provisioned IOPS."
  }
}

variable "bam_queue_enabled" {
  type        = bool
  default     = false
  description = "Enable separate BAM high (32-vCPU) and low (8-vCPU) fleets, one BAM task per worker."
}
variable "bam_max_vcpus" {
  type        = number
  default     = 256
  description = "BAM high capacity; retain the original variable name for existing tfvars."
  validation {
    condition     = var.bam_max_vcpus >= 32 && floor(var.bam_max_vcpus / 32) == var.bam_max_vcpus / 32
    error_message = "Use a positive multiple of 32 vCPUs for r8i.8xlarge BAM workers."
  }
}
variable "bam_low_max_vcpus" {
  type        = number
  default     = 64
  description = "BAM low capacity in vCPUs; each exclusive r8i.2xlarge worker reserves 8."
  validation {
    condition     = var.bam_low_max_vcpus >= 8 && floor(var.bam_low_max_vcpus / 8) == var.bam_low_max_vcpus / 8
    error_message = "BAM low capacity must be at least 8 and a multiple of 8 vCPUs."
  }
}
variable "bam_low_scratch_gb" {
  type        = number
  default     = 2000
  description = "BAM low gp3 root/staging disk in GiB; independent of high and CPU/GPU disks."
  validation {
    condition     = var.bam_low_scratch_gb >= 500 && var.bam_low_scratch_gb <= 16384 && floor(var.bam_low_scratch_gb) == var.bam_low_scratch_gb
    error_message = "BAM low scratch disk must be an integer from 500 to 16384 GiB."
  }
}
variable "bam_scratch_gb" {
  type        = number
  default     = 2000
  description = "BAM high gp3 root/staging disk in GiB; retain the original variable name for existing tfvars."
  validation {
    condition     = var.bam_scratch_gb >= 500 && var.bam_scratch_gb <= 16384 && floor(var.bam_scratch_gb) == var.bam_scratch_gb
    error_message = "Use an integer 500–16384 GiB."
  }
}

variable "cpu_ami_id" {
  type        = string
  description = "Pinned x86-64 ECS AL2023 AMI for CPU and BAM workers."
  validation {
    condition     = can(regex("^ami-[a-f0-9]{17}$", var.cpu_ami_id))
    error_message = "Supply an explicit ECS AL2023 AMI ID in the deployment region."
  }
}
variable "gpu_ami_id" {
  type        = string
  default     = null
  description = "Pinned x86-64 ECS AL2023 NVIDIA AMI for GPU workers."
  validation {
    condition     = (length(var.gpu_instance_types) == 0 && length(var.basecall_instance_types) == 0) || can(regex("^ami-[a-f0-9]{17}$", var.gpu_ami_id))
    error_message = "Supply an explicit NVIDIA AMI ID when GPU workers are enabled."
  }
}

variable "basecall_instance_types" {
  type        = list(string)
  default     = []
  description = "Dedicated x86-64 NVIDIA instance families for Dorado basecalling."
}
variable "basecall_max_vcpus" {
  type        = number
  default     = 32
  description = "Independent maximum capacity for the basecalling GPU fleet."
  validation {
    condition     = var.basecall_max_vcpus > 0
    error_message = "Basecalling GPU capacity must be positive."
  }
}
