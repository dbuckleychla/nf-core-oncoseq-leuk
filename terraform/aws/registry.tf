# Read the existing secret only on workers, never via a Terraform data source.
# The task and coordinator roles do not receive access to registry credentials.
resource "aws_iam_role_policy" "registry_auth" {
  name = "${var.name}-ghcr-pull"
  role = aws_iam_role.instance.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = concat([
      {
        Effect   = "Allow"
        Action   = ["secretsmanager:GetSecretValue"]
        Resource = [var.ghcr_secret_arn]
      }
      ], var.ghcr_secret_kms_key_arn == null ? [] : [
      {
        Effect   = "Allow"
        Action   = ["kms:Decrypt"]
        Resource = [var.ghcr_secret_kms_key_arn]
        Condition = {
          StringEquals = {
            "kms:ViaService"                  = "secretsmanager.${var.region}.${data.aws_partition.current.dns_suffix}"
            "kms:EncryptionContext:SecretARN" = var.ghcr_secret_arn
          }
        }
      }
    ])
  })
}
