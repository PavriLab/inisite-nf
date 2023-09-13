process CLUSTERSCAN_CLUSTER_IS {
    tag "meta.id"

    input:
    tuple val(meta), path(commonPeaks)

    output:
    tuple val(meta), path(commonPeaks),                 emit: peaks
    tuple val(meta), path("${meta.id}_clusters.bed"),   emit: clusters

    script:
    clusterScanScript = "${NXF_HOME}/assets/pavrilab/inisite-nf/bin/clusterscan.py"

    """
    clusterinitsites.py \
        -p !{commonPeaks} \
        --clusterScan !{clusterScanScript} \
        -o !{meta.id}
    """
}
