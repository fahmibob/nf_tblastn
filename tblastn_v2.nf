// enable dsl2
nextflow.enable.dsl=2
//specify output in parameter job and common option
params.job="default"
params.forking=7
params.filterbyFile=false

//tdecode
params.tdecode=false
params.dbName="$baseDir/tdecodeDB/ens_min_prot"

//makeblastdb option
params.tblastn=false
params.blastout=false
params.dustmask=false
params.windowmask=false

//tblastn option
params.query="$baseDir/query/compiledquery.fasta"
params.threads=4
params.outfmt=6

//count tblastn option
params.reference="$baseDir/ref/reftblastn.json"
params.module="$baseDir/bin/countBM.rb"
params.bitscore=50
params.countenzyme=true

workflow {
  if (params.tblastn) {
    Channel.fromPath(params.tblastn)
      .map {file -> tuple(file.simpleName, file)}
      .set{dbCH}

    if (params.tdecode){
      transdecoder(dbCH)
      dbCH=transdecoder.out
    }


    if (params.filterbyFile) {
      def filterList = new File(params.filterbyFile).collect {it}
      dbCH=dbCH.map { if (it[0].toString() in filterList){it}}
    }

    if (params.tblastn && params.dustmask) {
      unziping(dbCH)
      //uniq_seq(unziping.out)
      rundustmasker(unziping.out)
      tblastn(rundustmasker.out)
    } else if (params.tblastn && params.windowmask) {
      window_masking(dbCH)
      tblastn(window_masking.out)
    } else if (params.tblastn) {
      makeblastdbonly(dbCH)
      tblastn(makeblastdbonly.out)
    }

    if (params.countenzyme) {
      tblastnOut_count(tblastn.out)
      tblastnOut_count.out
        .collectFile(name:"${params.job}.txt", storeDir:"$baseDir/output_count/") {it}
    }
  } else if (params.blastout) {

    Channel.fromPath(params.tblastn)
      .map {file -> tuple(file.simpleName, file)}
      .set{blastoutCH}
    if (params.filterbyFile) {
      def filterList = new File(params.filterbyFile).collect {it}
      blastoutCH=blastoutCH.map { if (it[0].toString() in filterList){it}}
    }
    tblastnOut_count(blastoutCH)
    tblastnOut_count.out
      .collectFile(name:"${params.job}.txt", storeDir:"$baseDir/output_count/") {it}
  }

}

process transdecoder {
  maxForks params.forking
  errorStrategy 'ignore'
  input:
  tuple val(sampleName), path(data)

  output:
  tuple val(sampleName), path('*.transdecoder.cds')

  script:
  fullName=data.getName()
  """
  TransDecoder.LongOrfs -t $data
  blastp -query $fullName'.transdecoder_dir/longest_orfs.pep'\
  -db $params.dbName -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads $params.threads > $sampleName'.blastp'
  TransDecoder.Predict -t $data --retain_blastp_hits $sampleName'.blastp'
  """
}


process unziping {
  input:
  tuple val(sampleName), path(data)

  output:
  tuple val(sampleName), path('*ch.fasta')

  script:
  """
  if [[ $data == *.gz ]]
  then
    gunzip -c $data > $sampleName'ch.fasta'
  else
    cat $data > $sampleName'ch.fasta'
  fi
  """
}


process uniq_seq {
  input:
  tuple val(sampleName), path(data)

  output:
  tuple val(sampleName), path("*_U.fasta")

  script:
  """
  gawk '{if (\$0 ~/^>/) {h[\$1]++; \$1=\$1 "_" h[\$1]} print}' $data > $sampleName'_U.fasta'
  """
}

process window_masking {
  publishDir "$baseDir/DB/$params.job", mode: 'copy'
  maxForks params.forking
  errorStrategy 'ignore'

  input:
  tuple val(sampleName), path(data)

  output:
  tuple val(sampleName), path("*.n??")

  script:
  if (data.getExtension() == "gz"){
    checkedData=data.getBaseName
  } else {
    checkedData=data
  }
  """
  if [[ $data == *.gz ]]
  then
    gunzip $data
  fi

  windowmasker -in $checkedData -infmt fasta -mk_counts\
  -parse_seqids -out $sampleName'.counts'

  windowmasker -in $checkedData -infmt fasta -ustat $sampleName'.counts'\
  -outfmt maskinfo_asn1_bin -parse_seqids -out $sampleName'_mask.asnb'

  makeblastdb -in $checkedData -dbtype nucl -parse_seqids -mask_data $sampleName'_mask.asnb'\
  -out $sampleName -title $sampleName
  """
}

process rundustmasker {
  publishDir "$baseDir/DB/$params.job", mode: 'copy'
  
  maxForks params.forking

  input:
  tuple val(sampleName), path(data)

  output:
  tuple val(sampleName), path("*.n??")

  script:
  """
  dustmasker -in $data -infmt fasta -parse_seqids \
  -outfmt maskinfo_asn1_bin -out $sampleName'_mask.asnb'

  makeblastdb -in $data -dbtype nucl -parse_seqids -mask_data $sampleName'_mask.asnb'\
  -out $sampleName -title $sampleName
   """
}

process masking_with_lowercase {
  publishDir "$baseDir/DB/$params.job", mode: 'copy'
  errorStrategy 'ignore'
  maxForks params.forking

  input:
  tuple val(sampleName), path(data)

  output:
  tuple val(sampleName), path("*.n??")

  script:
  """
  if [[ $data == *.gz ]]
  then
    gunzip -c $data | convert2blastmask -in - -parse_seqids -masking_algorithm repeat\
    -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out $sampleName'_mask.asnb'

    gunzip -c $data | makeblastdb -in - -dbtype nucl -parse_seqids -mask_data $sampleName'_mask.asnb'\
    -out $sampleName -title $sampleName
  else
    convert2blastmask -in $data -parse_seqids -masking_algorithm repeat\
    -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out $sampleName'_mask.asnb'

    makeblastdb -in $data -dbtype nucl -parse_seqids -mask_data $sampleName'_mask.asnb'\
    -out $sampleName -title $sampleName
  fi
  """
}

process makeblastdbonly {
  publishDir "$baseDir/DB/$params.job", mode: 'copy'
  maxForks params.forking

  input:
  tuple val(sampleName), path(data)

  output:
  tuple val(sampleName), path("*.n??")

  script:
  """
  if [[ $data == *.gz ]]
  then
    gunzip -c $data | makeblastdb -in - -dbtype nucl -parse_seqids\
    -out $sampleName -title $sampleName
  else
    makeblastdb -in $data -dbtype nucl -parse_seqids\
    -out $sampleName -title $sampleName
  fi
  """

}

process tblastn {
  publishDir "$baseDir/output_blast/$params.job", mode: 'copy'
  errorStrategy 'ignore'
  maxForks params.forking

  input:
  tuple val(dbName), path(data)

  output:
  tuple val(dbName), path("*.blastout")

  script:
  """
  tblastn -num_threads $params.threads -query $params.query -db $dbName\
  -out $dbName'.blastout' -outfmt $params.outfmt
  """
}


process tblastnOut_count {
  errorStrategy 'ignore'
  maxForks params.forking
  input:
  tuple val(sampleName), path(data)

  output:
  path("*.out")

  script:
  """
  ruby -r '$params.module' -e 'filter_identity("$data", "$params.reference", $params.bitscore, "tblastn")' > $sampleName'.out'
  """
}
