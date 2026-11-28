process CRAMINO {
    output:
    stdout emit: settings
    script:
    def exitCode = task.attempt == 1 ? params.test_exit_code : 0
    """
    printf '%s\\n' '${task.queue}' '${task.cpus}' '${task.memory}' '${task.attempt}'
    exit ${exitCode}
    """
}

process SAMTOOLS_COUNT_READS {
    output:
    stdout emit: settings
    script:
    """
    printf '%s\\n' '${task.attempt}'
    exit 137
    """
}
