
rule associated_genes_from_genebass:
    conda: '../../envs/ukbb-splicing-analyisi-hail.yaml'
    input: genebass = '/s/project/bayesRare/rare_common_analysis/genebass/results.mt',
           traits = config['dataDir'] + '/input/traits_monti.tsv'
    output: filtered_genebass = config['dataDir'] + '/genebass/filtered_associations.tsv'
    script: '../scripts/extractBurdenResults.py'