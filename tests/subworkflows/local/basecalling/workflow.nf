include { BASECALL_SHARDED } from '../../../../subworkflows/local/basecalling/basecall_sharded'

workflow CHECK_SHARDS {

    main:
    def repo = params.shard_test_root ?: projectDir
    def root = file("${repo}/tests/subworkflows/local/basecalling/fixtures")
    def model = root.resolve('base_model')
    def empty = file("${repo}/assets/NOFILE")
    def mod = file("${repo}/assets/NOMOD")
    def rows = params.single_pod5
        ? [[[id: 'alpha'], root.resolve('alpha/a.pod5'), empty, model, mod]]
        : [[[id: 'alpha'], root.resolve('alpha'), empty, model, mod], [[id: 'beta'], root.resolve('beta'), empty, model, mod]]
    BASECALL_SHARDED(channel.fromList(rows))

    emit:
    ubam = BASECALL_SHARDED.out.ubam
}
