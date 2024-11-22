process concat_csv {
    container "${params.container__pandas}"
    label 'io_limited'
    publishDir "${params.output}/${params.output_subfolder}/", mode: 'copy', overwrite: true, enabled: params.publish
    
    input:
    path "*.csv"
    
    output:
    path "${params.output_csv}"
    
    script:
    template "concat_csv.py"
}