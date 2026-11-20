process DORADO_BASECALL {
    label 'process_gpu'
    output:
    stdout emit: settings
    script:
    """
    printf '%s\\n' '${task.queue}' '${task.accelerator.request}' '${task.container}' '${task.containerOptions}' '${task.clusterOptions}'
    """
}

process CLASSY_COMBINED {
    output:
    stdout emit: settings
    script:
    """
    printf '%s\\n' '${task.queue}' '${task.cpus}' '${task.container}'
    """
}

process BCFTOOLS_FILTER {
    output:
    stdout emit: settings
    script:
    """
    printf '%s\\n' '${task.container}'
    """
}
