#!/usr/bin/env nextflow

if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 
if (!params.adapter){params.adapter = ""}
read_ch = file(params.adapter, type: 'any')

if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Nem található egyező leolvasás: ${params.reads}" }
	.into{reads_in_1;reads_in_2}
 } else {  
	reads_in_1 = Channel.empty()
	reads_in_2 = Channel.empty()
 }

Channel.value(params.mate).into{mate_1;mate_2;mate_3;mate_4;mate_5;mate_6;mate_7}

process Adapter_Removal {

		publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /reads\/.*.fastq$/) "$name/adapter_removal/$filename"}
		input:
		 set val(name), file(reads) from reads_in_2
		 val mate from mate_4

		output:
		 set val(name), file("reads/*.fastq")  into reads_out_1
		 file "*.{fastx,trimmomatic}.log"  into log_file_1

		errorStrategy 'retry'

when:
		(params.run_Adapter_Removal && (params.run_Adapter_Removal == "yes")) || !params.run_Adapter_Removal

shell:
	phred = params.Adapter_Removal.phred
	Tool_for_Adapter_Removal = params.Adapter_Removal.Tool_for_Adapter_Removal
	Adapter_Sequence = params.Adapter_Removal.Adapter_Sequence
//trimmomatic_inputs
	min_length = params.Adapter_Removal.min_length
	seed_mismatches = params.Adapter_Removal.seed_mismatches
	palindrome_clip_threshold = params.Adapter_Removal.palindrome_clip_threshold
	simple_clip_threshold = params.Adapter_Removal.simple_clip_threshold

//fastx_clipper_inputs
	discard_non_clipped = params.Adapter_Removal.discard_non_clipped
    
remove_previous_reads = params.Adapter_Removal.remove_previous_reads
discard_non_clipped_text = ""
if (discard_non_clipped == "yes") {discard_non_clipped_text = "-c"}
nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage;
 use Cwd qw();
 
runCmd("mkdir reads adapter unpaired");

open(OUT, ">adapter/adapter.fa");
my @adaps=split(/\n/,"!{Adapter_Sequence}");
my $i=1;
foreach my $adap (@adaps)
{
 print OUT ">adapter$i\\n$adap\\n";
 $i++;
}
close(OUT);

system("!{runGzip}");
my $quality="!{phred}";
print "fastq quality: $quality\\n";
print "tool: !{Tool_for_Adapter_Removal}\\n";

if ("!{mate}" eq "pair") {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic PE -threads 1 -phred${quality} !{file1} !{file2} reads/!{name}.1.fastq unpaired/!{name}.1.fastq.unpaired reads/!{name}.2.fastq unpaired/!{name}.2.fastq.unpaired ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        print "Fastx_clipper használata nem javasolt pair fájlok esetén.";
    }
} else {
    if ("!{Tool_for_Adapter_Removal}" eq "trimmomatic") {
        runCmd("trimmomatic SE -threads 1  -phred${quality} !{file1} reads/!{name}.fastq ILLUMINACLIP:adapter/adapter.fa:!{seed_mismatches}:!{palindrome_clip_threshold}:!{simple_clip_threshold} MINLEN:!{min_length} 2> !{name}.trimmomatic.log");
    } elsif ("!{Tool_for_Adapter_Removal}" eq "fastx_clipper") {
        runCmd("fastx_clipper  -Q $quality -a !{Adapter_Sequence} -l !{min_length} !{discard_non_clipped_text} -v -i !{file1} -o reads/!{name}.fastq > !{name}.fastx.log");
    }
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir = join '/', @paths;
    $inputsdir .= "/work";
    print "INFO: Leolvasási bemenetek törlődnek, ha a $workdir $inputsdir ba találhatóak\\n";
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile törölve.\\n";
        }
    }
}


##Subroutines
sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}
'''

}

process Trimmer {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /reads\/.*q$/) "$name/trimmer/$filename"}
input:
 set val(name), file(reads) from reads_out_1
 val mate from mate_5

output:
 set val(name), file("reads/*q")  into reads_out_trimmmer
 file "*.log" optional true  into log_file_trimmer

errorStrategy 'retry'

when:
(params.run_Trimmer && (params.run_Trimmer == "yes")) || !params.run_Trimmer

shell:
phred = params.Trimmer.phred
single_or_paired_end_reads = params.Trimmer.single_or_paired_end_reads
trim_length_5prime = params.Trimmer.trim_length_5prime
trim_length_3prime = params.Trimmer.trim_length_3prime
trim_length_5prime_R1 = params.Trimmer.trim_length_5prime_R1
trim_length_3prime_R1 = params.Trimmer.trim_length_3prime_R1
trim_length_5prime_R2 = params.Trimmer.trim_length_5prime_R2
trim_length_3prime_R2 = params.Trimmer.trim_length_3prime_R2
remove_previous_reads = params.Trimmer.remove_previous_reads

nameAll = reads.toString()
nameArray = nameAll.split(' ')
file2 = ""
if (nameAll.contains('.gz')) {
    newName =  nameArray[0] - ~/(\.fastq.gz)?(\.fq.gz)?$/
    file1 =  nameArray[0] - '.gz' 
    if (mate == "pair") {file2 =  nameArray[1] - '.gz'}
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    newName =  nameArray[0] - ~/(\.fastq)?(\.fq)?$/
    file1 =  nameArray[0]
    if (mate == "pair") {file2 =  nameArray[1]}
    runGzip = ''
}
'''
#!/usr/bin/env perl
 use List::Util qw[min max];
 use strict;
 use File::Basename;
 use Getopt::Long;
 use Pod::Usage; 
 use Cwd qw();
 
system("mkdir reads");
system("!{runGzip}");
my $file1 = "";
my $file2 = "";
if ("!{mate}" eq "pair") {
    $file1 = "!{file1}";
    $file2 = "!{file2}";
    my $trim1 = "!{trim_length_5prime_R1}:!{trim_length_3prime_R1}";
    my $trim2 = "!{trim_length_5prime_R2}:!{trim_length_3prime_R2}";
    my $len=getLength($file1);
    print "length of $file1: $len\\n";
    trimFiles($file1, $trim1, $len);
    my $len=getLength($file2);
    print "INFO: length of $file2: $len\\n";
    trimFiles($file2, $trim2, $len);
} else {
    $file1 = "!{file1}";
    my $trim1 = "!{trim_length_5prime}:!{trim_length_3prime}";
    my $len=getLength($file1);
    print "INFO: length of file1: $len\\n";
    trimFiles($file1, $trim1, $len);
}
if ("!{remove_previous_reads}" eq "true") {
    my $currpath = Cwd::cwd();
    my @paths = (split '/', $currpath);
    splice(@paths, -2);
    my $workdir= join '/', @paths;
    splice(@paths, -1);
    my $inputsdir= join '/', @paths;
    $inputsdir .= "/inputs";
    print "INFO: Leolvasási bemenetek törlődnek, ha a workdir inputsdir ba találhatóak\\n"; 
    my @listOfFiles = `readlink -e !{file1} !{file2}`;
    foreach my $targetFile (@listOfFiles){
        if (index($targetFile, $workdir) != -1 || index($targetFile, $inputsdir) != -1) {
            system("rm -f $targetFile");
            print "INFO: $targetFile törölve.\\n";
        }
    }
}



sub trimFiles
{
  my ($file, $trim, $len)=@_;
    my @nts=split(/[,:\\s\\t]+/,$trim);
    my $inpfile="";
    my $com="";
    my $i=1;
    my $outfile="";
    my $param="";
    my $quality="-Q!{phred}";

    if (scalar(@nts)==2)
    {
      $param = "-f ".($nts[0]+1) if (exists($nts[0]) && $nts[0] >= 0 );
      $param .= " -l ".($len-$nts[1]) if (exists($nts[0]) && $nts[1] > 0 );
      $outfile="reads/$file";  
      $com="fastx_trimmer $quality -v $param -o $outfile -i $file > !{name}.fastx_trimmer.log" if ((exists($nts[0]) && $nts[0] > 0) || (exists($nts[0]) && $nts[1] > 0 ));
      print "INFO: $com\\n";
      if ($com eq ""){
          print "INFO: Trimmer átugrotta a $file fájlt.\\n";
          system("mv $file reads/.");
      } else {
          runCmd("$com");
          print "INFO: Trimmer elvégzi a $file fájlt elemzését.\\n";
      }
    }

    
}


sub getLength
{
   my ($filename)=@_;
   open (IN, $filename);
   my $j=1;
   my $len=0;
   while(my $line=<IN>)
   {
     chomp($line);
     if ($j >50) { last;}
     if ($j%4==0)
     {
        $len=length($line);
     }
     $j++;
   }
   close(IN);
   return $len;
}

sub runCmd {
    my ($com) = @_;
    my $error = system($com);
    if   ($error) { die "Command failed: $error $com\\n"; }
    else          { print "Command successful: $com\\n"; }
}

'''

}

process Trimmer_Removal_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /$name.trimmer_summary.tsv$/) "$name/trimmer_summary/$filename"}
input:
 file logfile from log_file_trimmer.collect()
 val mate from mate_7

output:
 file "$name.trimmer_summary.tsv"  into output_tsv_rem

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx_trimmer\\.log/){
        $file =~ /(.*)\\.fastx_trimmer\\.log/;
        my $mapper   = "fastx_trimmer";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        my $in;
        my $out;
        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Trimmer" ];
    }
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "$name.trimmer_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

process Adapter_Removal_Summary {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /adapter_removal_summary.tsv$/) "$name/adapter_removal_summary/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /adapter_removal_detailed_summary.tsv$/) "$name/adapter_removal_detailed_summary/$filename"}
input:
 file logfile from log_file_1.collect()
 val mate from mate_2

output:
 file "adapter_removal_summary.tsv"  into output_tsv_rem_1
 file "adapter_removal_detailed_summary.tsv" optional true  into output_tsv_det

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %tsvDetail;
my %headerHash;
my %headerText;
my %headerTextDetail;

my $i = 0;
chomp( my $contents = `ls *.log` );

my @files = split( /[\\n]+/, $contents );
foreach my $file (@files) {
    $i++;
    my $mapOrder = "1";
    if ($file =~ /(.*)\\.fastx\\.log/){
        $file =~ /(.*)\\.fastx\\.log/;
        my $mapper   = "fastx";
        my $name = $1;    ##sample name
        push( @header, $mapper );

        my $in;
        my $out;
        my $tooshort;
        my $adapteronly;
        my $noncliped;
        my $Nreads;

        chomp( $in =`cat $file | grep 'Input:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $out =`cat $file | grep 'Output:' | awk '{sum+=\\$2} END {print sum}'` );
        chomp( $tooshort =`cat $file | grep 'too-short reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $adapteronly =`cat $file | grep 'adapter-only reads' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $noncliped =`cat $file | grep 'non-clipped reads.' | awk '{sum+=\\$2} END {print sum}'`);
        chomp( $Nreads =`cat $file | grep 'N reads.' | awk '{sum+=\\$2} END {print sum}'` );

        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        $tsvDetail{$name}{$mapper} = [ $in, $tooshort, $adapteronly, $noncliped, $Nreads, $out ];
        $headerTextDetail{$mapOrder} = ["Total Reads","Too-short reads","Adapter-only reads","Non-clipped reads","N reads","Reads After Adapter Removal"];
    } elsif ($file =~ /(.*)\\.trimmomatic\\.log/){
        $file =~ /(.*)\\.trimmomatic\\.log/;
        my $mapper   = "trimmomatic";
        my $name = $1;    ##sample name
        push( @header, $mapper );
        
        my $in;
        my $out;

        if ( "!{mate}" eq "pair"){
            chomp( $in =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$4} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Read Pairs:' | awk '{sum+=\\$7} END {print sum}'` );
        } else {
            chomp( $in =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$3} END {print sum}'` );
            chomp( $out =`cat $file | grep 'Input Reads:' | awk '{sum+=\\$5} END {print sum}'` );
        }
        


        $tsv{$name}{$mapper} = [ $in, $out ];
        $headerHash{$mapOrder} = $mapper;
        $headerText{$mapOrder} = [ "Total Reads", "Reads After Adapter Removal" ];
        
    }
    
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary          = "adapter_removal_summary.tsv";
my $detailed_summary = "adapter_removal_detailed_summary.tsv";
writeFile( $summary,          \\%headerText,       \\%tsv );
if (%headerTextDetail){
    writeFile( $detailed_summary, \\%headerTextDetail, \\%tsvDetail );  
}

sub writeFile {
    my $summary    = $_[0];
    my %headerText = %{ $_[1] };
    my %tsv        = %{ $_[2] };
    open( OUT, ">$summary" );
    print OUT "Sample\\t";
    my @headArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push( @headArr, @{ $headerText{$mapOrder} } );
    }
    my $headArrAll = join( "\\t", @headArr );
    print OUT "$headArrAll\\n";

    foreach my $name ( keys %tsv ) {
        my @rowArr = ();
        for my $mapOrder (@sortedOrderArray) {
            push( @rowArr, @{ $tsv{$name}{ $headerHash{$mapOrder} } } );
        }
        my $rowArrAll = join( "\\t", @rowArr );
        print OUT "$name\\t$rowArrAll\\n";
    }
    close(OUT);
}

'''
}

process FastQC {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.(html|zip)$/) "$name/fastqc/$filename"}
input:
 val mate from mate_1
 set val(name), file(reads) from reads_in_1

output:
 file '*.{html,zip}'  into fastq_out

errorStrategy 'retry'
maxRetries 3

when:
(params.run_FastQC && (params.run_FastQC == "yes"))

script:
nameAll = reads.toString()
if (nameAll.contains('.gz')) {
    file =  nameAll - '.gz' - '.gz'
    runGzip = "ls *.gz | xargs -i echo gzip -df {} | sh"
} else {
    file =  nameAll 
    runGzip = ''
}
"""
${runGzip}
fastqc ${file} 
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}