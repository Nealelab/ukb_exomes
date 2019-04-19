import hail as hl

bucket = 'gs://ukbb-pharma-exome-analysis'
root = f'{bucket}/mt'
# common_high_callrate_path = f'{root}/ukb.common.high.callrate.mt'
# common_high_callrate_pruned_path = f'{root}/ukb.common.high.callrate.pruned.mt'
common_high_callrate_pruned_plink_path = f'{root}/ukb.common.high.callrate.pruned.plink'

public_bucket = 'gs://ukbb-exome-public/'
