require 'json'
def ranges_overlap?(range_a, range_b)
  range_b.begin <= range_a.end && range_a.begin <= range_b.end
end

#filter overlapped contig based on its match position
def filter_overlapped_contigs(aryblasout,blast_type)
  grouped_contigs=[] #create sub array of match query
  aryblasout.chunk {|f|
    f.split("\t")[0]
  }.each {|ch, group| grouped_contigs.push(group)}

  filoverlapres=[]
  #filtering matched query with overlapped length
  grouped_contigs.map do |sub_ary|
    unless sub_ary.count < 2
      if blast_type == "blastx"
        sortedrange=sub_ary.map {|f|(f.split("\t")[6]).to_i..(f.split("\t")[7]).to_i}
      else
        sortedrange=sub_ary.map {|f|(f.split("\t")[8]).to_i..(f.split("\t")[9]).to_i}
      end
      sortedrange=sortedrange.map {|f| f.begin > f.end ? (f.end)..(f.begin) : (f.begin)..(f.end)}
      mod_sortedrange=sortedrange.uniq
      for i in 0..mod_sortedrange.count-1
        deleted_list=[]
        for j in i+1..mod_sortedrange.count-1
          if ranges_overlap?(mod_sortedrange[i], mod_sortedrange[j])
            deleted_list.push(mod_sortedrange[j])
          end
        end
        if deleted_list.any?
          deleted_list.each {|del|mod_sortedrange.delete(del)}
        end
      end
      mod_sortedrange.map!{|f|sortedrange.index(f)}
      mod_sortedrange.each {|f| filoverlapres.push(sub_ary[f])}
    else
      filoverlapres.push(sub_ary[0])
    end
  end
  return filoverlapres
end


def calculate_enzymes(inputData,reference=nil)
  unless reference == nil
    enzyme_ref=[]
    enzyme_count=[]
    for i in 0..reference.keys.count-1
      enzyme_ref.push("#{reference.keys[i]}")
      enzyme_count.push(0)
    end
    for i in 0..inputData.count-1
      unless (reference.find {|key, value| value.include?(inputData[i].split("\t")[0])}).nil?
        enzyme_count[enzyme_ref.index(reference.find {|key, value| value.include?(inputData[i].split("\t")[0])}.first)]+=1
      end
    end
    return enzyme_count
  else
    result=["#{inputData.count}"]
    return result
  end

end

#filter blastout using set of identity from json file and bitscore
def filter_identity(blastout, bitscore, blast_type, threshold_id=30, jsonfile=nil)
  filename=File.basename(blastout, ".blastout")
  unless jsonfile == nil
    reference = JSON.parse(File.read(jsonfile))
    idFilteredData=File.readlines(blastout).select {|line|
      unless line[0] == '#'
        if blast_type == "blastx" && reference.find { |key, values| values.include?(line.split("\t")[1])} != nil
          line.split("\t")[2].to_f > (reference.find { |key, values| values.include?(line.split("\t")[1])}.first).split("_")[1].to_f && line.split("\t")[11].to_i > bitscore
        elsif reference.find { |key, values| values.include?(line.split("\t")[0])} != nil
          line.split("\t")[2].to_f > (reference.find { |key, values| values.include?(line.split("\t")[0])}.first).split("_")[1].to_f && line.split("\t")[11].to_i > bitscore
        end
      end
    }
    ovFilteredData=filter_overlapped_contigs(idFilteredData, blast_type)
    enzymeAmount=calculate_enzymes(ovFilteredData, reference)
  else
    idFilteredData=File.readlines(blastout).select {|line|
      unless line[0] == '#'
        line.split("\t")[2].to_f > threshold_id && line.split("\t")[11].to_i > bitscore
      end
    }
    ovFilteredData=filter_overlapped_contigs(idFilteredData, blast_type)
    enzymeAmount=calculate_enzymes(ovFilteredData)
  end
  puts "#{filename}\t#{enzymeAmount.join("\t")}\t#{ovFilteredData.count}"
end
