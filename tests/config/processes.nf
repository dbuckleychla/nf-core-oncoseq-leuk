process DORADO_BASECALL {
    label 'process_gpu'
    output:
    stdout emit: settings
    script:
    """
    printf '%s\\n' '${task.queue}' '${task.accelerator.request}' '${task.container}' '${task.containerOptions}' '${task.clusterOptions}' '${task.cpus}' '${task.memory}' '${task.maxRetries}'
    """
}

process CLASSY_COMBINED {
    output:
    stdout emit: settings
    script:
    """
    printf '%s\\n' '${task.queue}' '${task.cpus}' '${task.container}' '${task.memory}' '${task.time}' '${task.ext.threads}'
    """
}

process SAMTOOLS_VIEW {
    output:
    stdout emit: settings
    script:
    """
    printf '%s\\n' '${task.queue}' '${task.cpus}' '${task.container}'
    """
}

process ENSEMBLVEP_VEP {
    output:
    stdout emit: settings
    script:
    """
    printf '%s\\n' '${task.queue}' '${task.cpus}' '${task.memory}' '${task.time}'
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
