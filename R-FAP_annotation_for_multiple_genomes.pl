#!/usr/bin/perl

$start_time = localtime(time);
print "Started at $start_time\n\n";

use Getopt::Long;
use lib 'lib';
use String::Simrank;
use autodie;
use File::Path qw(make_path remove_tree);

print "\n--------------------------------------------------------------------------------\n";
print "Welcome to the main Rapid Functional Annotation of Prokaryotes (R-FAP) tool.\n";
print "database builder. This software will take the complete gene translations from a\n";
print "set of selected strains and provide annotations for them in under a minute.\n";
print "--------------------------------------------------------------------------------\n\n\n";

print "Please provide the full path to the directory of the .faa files to be annotated:\n";

$operating_directory = <STDIN>;
chomp $operationg_directory;

if($operating_directory !~ /\w/) #exit if user provides no directory
{	
    	die "You must specify where to find the input files!\n\n"; 
}
elsif($operating_directory !~ /^\// || $operating_directory !~ /\/$/)
{#in case user forgets to include the '/' at the beginning or end
	if($operating_directory !~ /^\//)
	{		
		substr($operating_directory,0,0) = "/";
	}

	if($operating_directory !~ /\/$/)
	{
		$operating_directory .= "/";
	}
}
chomp $operating_directory;

#these two 'make_path' commands create the subdirectories for the output files
make_path($operating_directory."new_annotation_files/");
make_path($operating_directory."new_unmatched_sequence_files");

print "Next, please enter the full path to the FASTA-formatted pan-genome database\n";
print "which will be used to match query proteins to their appropriate annotations:\n\n";

$pangenome_filepath = <STDIN>;

if($pangenome_filepath !~ /\w/) 
{	
    	die "You must specify where to find the database files!\n\n"; 
}
elsif($pangenome_filepath !~ /^\//)
{		
	substr($pangenome_filepath,0,0) = "/";
}

chomp $pangenome_filepath;

#if($pangenome_filepath =~ /\n\z/)
#{
#    chop($pangenome_filepath);
#}
print "\nNow, please specify the full path to the three-column annotation source\n";
print "file. (It would probably be best to make sure its located in the operating\n";
print "directory, but in either case you need to enter the full path to it here.)\n\n";

$annotation_filename = <STDIN>;

if($annotation_filename !~ /\w/)
{
    die "You must provide a source of annotations!\n\n";
}
elsif($annotation_filename !~ /^\//)
{		
    substr($annotation_filename,0,0) = "/";
}

open(ANNO,$annotation_filename) || die "Could not open $annotation_filename\n\n";
$anno_iterator = 0;

# This loop iterates through the annotation file, which is a three-column, tab-delimited
# file of basic annotation for each gene and the tool used to determine it
while($buf = <ANNO>)
{
    chomp $buf;
    if($anno_iterator > 0)
    {
	($protein,$id_tool,$annotation) = split(/\t/, $buf);
	
	if($id_tool !~ /\w/)
	{
	    $id_tool = "";
	}

	$Annotation_hash{$protein} = "$id_tool\t$annotation";
    }
    $anno_iterator++;
}

# Below are the commands which declare a new Simrank object and create a new binary database file
$sr = new String::Simrank ({ data => $pangenome_filepath });
if ($cl_args->{"rebuild"} || !$sr->{binary_ready} ) 
{
    $sr->formatdb({ wordlen => 10,
		    minlen  => 10,
		  });
}

$i=0;


# This loop will iterate through all of the '.faa' files in the specified directory and annotate them
@files = <$operating_directory*.faa>;
$number_of_files = scalar(@files);
print "\nThe number of .faa files in $operating_directory is $number_of_files\n";

#while (defined($sequence_filename = glob '"$operating_directory"."*.faa"' )) 
foreach $sequence_filename (@files)
{
   
    open(SEQ,$sequence_filename) || die "Did not open $sequence_filename.\n\n";
        
    $sequence_name = '';
    
    while($buf = <SEQ>)
    {
	chomp $buf;
	
	if($buf =~ /\>/)
	{	    
	    $sequence_name = substr $buf, 1;
	    $sequence = '';
	    $length = 0;
	}
	elsif($buf !~ /\>/ && $buf =~ /\w/)
	{
	    $Sequences{$sequence_name} .= $buf;
	}
	else
	{
	}
    }
    close SEQ;
    
    $number_of_query_genes = scalar(keys %Sequences);
    $current_time = localtime(time);
    print "\nStarted at $start_time\nCurrent time is $current_time\n\n";
    print "Number of genes in $sequence_filename = $number_of_query_genes\n";
    
    #30% appears to be the lowest point at which accuracy is still maintained at 100%
    $minimum_match_percentage = 30;

    #this is actually just a count of the number of secondary matches found by Simrank
    $number_of_false_positives = 0; 

    #this is actually just a percentage of the total proteins in a genome which Simrank doesn't find
    $false_negative_rate = 0; 
    
    $j=0;
    $number_of_false_positives = 0;
    undef %Match_list;
    undef @putative_new_genes;
    
    # The following command calculates the match between each query protein and each sequence
    # in the search database; these matches are sorted most- to least-similar and retunred as
    # a hash (specifically, a hash of arrays)
    $matches = $sr->match_oligos( { query => $sequence_filename,
				    outlen => 2,
				    minpct => $minimum_match_percentage,
				    silent => true,
				    valid_chars => 'ABCDEFGHIJKLMNOPQRSTUVWYXZ',
				  });
    
    # Once this hash has been created, one must step through it protein at a time; each query
    # protein has an entry in the hash, even if it didn't match anything in the database.
    foreach $key (keys %{$matches})
    {
	#this 'if' tests to see if a protein's entry in the hash is blank; if it is, then it
	#didn't match to anything in the database above the specified similarity threshold.
	if( @{ $matches->{$key} } == '')
	{
	    push(@putative_new_genes, $key);
	}
	
	$k = 0;
	foreach $hit ( @{ $matches->{$key} } ) 
	{
	    $gene_name = $key;
	    $match = $hit->[0];
	    $match_percent = $hit->[1];
	    
	    if($k==0) #this indicates the first (and therefore best) match
	    {
		$Match_list{$gene_name} .= "$gene_name\t$Annotation_hash{$match}";
	    }
	    elsif($k>0) #even if this is ever true, it represents a suprious match
	    {
		$number_of_false_positives++
	    }
	    
	    $k++;
	}
	$j++;
    }
    
    $number_of_matched_genes = scalar(keys %Match_list);
    
    #the following calculates the 'miss rate' by dividing the number which ACTUALLY matched
    #at or above the given threshold from the total number of proteins searched.
    $false_negative_rate = 100*(1-($number_of_matched_genes/$number_of_query_genes));
    $i++;
    print "Currently on genome $i of $number_of_files\n\tMissed rate = $false_negative_rate\n";
    print "\tNumber of secondary matches = $number_of_false_positives\n";
    
    $short_input_filename_start = rindex($sequence_filename,"/");
    $short_input_filename = substr($sequence_filename,$short_input_filename_start);
    
    while($short_input_filename =~ /\./)
    {
	chop $short_input_filename;
    }
    
    #create the output file from the input file name
    $final_output_filename = $operating_directory."new_annotation_files".$short_input_filename."_annotation.txt";
    open(OUT, ">$final_output_filename") || die "Could not create and/or open $final_output_filename.\n\n";
    
    #create the name for the FASTA file of unmatched proteins
    $unmatched_genes_filename = $operating_directory."new_unmatched_sequence_files".$short_input_filename."_unmatched.fasta";
    open(UNMATCHED, ">$unmatched_genes_filename") || die "Could not create and/or open $unmatched_genes_filename.\n\n";
    
    print OUT "Protein\tID_tool\tAnnotation\n";
    foreach $gene (keys %Match_list)
    {
	print OUT "$Match_list{$gene}\n"; #note that this is formatted like the annotation file
    }
    
    close OUT;
    
    
    for($k=0;$k<scalar(@putative_new_genes);$k++)
    {
	print UNMATCHED ">$putative_new_genes[$k]\n$Sequences{$putative_new_genes[$k]}\n";
    }

    close UNMATCHED;
    
    #it's important to 'undef' (erase) these data structures each time
    undef @putative_new_genes;
    undef %Match_list;
    undef %Sequences;

}
    
$current_time = localtime(time);
print "\n\nStarted at $start_time\nFinished at $current_time\n\n\n";
    
