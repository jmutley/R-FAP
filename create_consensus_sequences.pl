#!/usr/bin/perl
$start_time = localtime;
print "$start_time\n\n";


print "The input file should contain one or more multiple sequence\nalignments (MSA) and must be in CLUSTAL (.clw) format\n\nPlease enter the full path to the file:\n\n";

$in = <STDIN>;
open(IN, $in) || die "Could not open $in.\n\n";

print "\n\nThe output file will be a FASTA-formatted file containing\nthe consensus sequence from each MSA in the input file,\ndetermined using the majority nucleotide or residue at each position.\n";
print "\nPlease specify the full path to the output file, whether new or existing:\n\n";

$out = <STDIN>;
open(OUT, ">$out") || die "Could not creat and/or open $out.\n\n";
$iterator = 0;
$gene_name = '';
$protein_sequence = '';
$consensus = '';
$consensus_sequence = '';
$gene_counter = 0;

$A = "A";
$B = "B";
$C = "C";
$D = "D";
$E = "E";
$F = "F";
$G = "G";
$H = "H";
$I = "I";
$J = "J";
$K = "K";
$L = "L";
$M = "M";
$N = "N";
$O = "O";
$P = "P";
$Q = "Q";
$R = "R";
$S = "S";
$T = "T";
$U = "U";
$V = "V";
$W = "W";
$X = "X";
$Y = "Y";
$Z = "Z";
$gap = "";


while($buf = <IN>)
{
	chomp ($buf);

	if($buf =~ /multiple/ || eof) #indicates a new alignment in the input file (or that the next line will be the end of the file)
	{
		if($iterator > 0)
		{
			$length = length($Protein_hash{$gene_name});
			for($i=0;$i<$length;$i++)
			{
				$consensus_character = substr($consensus, $i, 1);
				


				if($consensus_character eq "*") #indicates a completely conserved position
				{
				    $sequence_char = substr($Protein_hash{$gene_name}, $i, 1);
				    $consensus_sequence .= $sequence_char;
				}
				else #indicates a position which is not conserved
				{
				    undef $Char_hash;

				    $a_counter = 0;
				    $b_counter = 0;
				    $c_counter = 0;
				    $d_counter = 0;
				    $e_counter = 0;
				    $f_counter = 0;
				    $g_counter = 0;
				    $h_counter = 0;
				    $i_counter = 0;
				    $j_counter = 0;
				    $k_counter = 0;
				    $l_counter = 0;
				    $m_counter = 0;
				    $n_counter = 0;
				    $o_counter = 0;
				    $p_counter = 0;
				    $q_counter = 0;
				    $r_counter = 0;
				    $s_counter = 0;
				    $t_counter = 0;
				    $u_counter = 0;
				    $v_counter = 0;
				    $w_counter = 0;
				    $x_counter = 0;
				    $y_counter = 0;
				    $z_counter = 0;
				    $gap_counter = 0;
				    
				    #iterate through each protein or DNA sequence to determine the most common residue at the given position
				    foreach $protein (keys %Protein_hash) 
				    {
					$consensus_char = substr($Protein_hash{$protein}, $i, 1);
					
					#this 'if/elsif' counts the number of times each residue occurs at the given position
					if($consensus_char eq "A" || $consensus_char eq "a")
					{
					    $a_counter++;
					}
					elsif($consensus_char eq "B" || $consensus_char eq "b")
					{
					    $b_counter++;
					}
					elsif($consensus_char eq "C" || $consensus_char eq "c")
					{
					    $c_counter++;
					}
					elsif($consensus_char eq "D" || $consensus_char eq "d")
					{
					    $d_counter++;
					}
					elsif($consensus_char eq "E" || $consensus_char eq "e")
					{
					    $e_counter++;
					}
					elsif($consensus_char eq "F" || $consensus_char eq "f")
					{
					    $f_counter++;
					}
					elsif($consensus_char eq "G" || $consensus_char eq "g")
					{
					    $g_counter++;
					}
					elsif($consensus_char eq "H" || $consensus_char eq "h")
					{
					    $h_counter++;
					}
					elsif($consensus_char eq "I" || $consensus_char eq "i")
					{
					    $i_counter++;
					}
					elsif($consensus_char eq "J" || $consensus_char eq "j")
					{
					    $j_counter++;
					}
					elsif($consensus_char eq "K" || $consensus_char eq "k")
					{
					    $k_counter++;
					}
					elsif($consensus_char eq "L" || $consensus_char eq "l")
					{
					    $l_counter++;
					}
					elsif($consensus_char eq "M" || $consensus_char eq "m")
					{
					    $m_counter++;
					}
					elsif($consensus_char eq "N" || $consensus_char eq "n")
					{
					    $n_counter++;
					}
					elsif($consensus_char eq "O" || $consensus_char eq "o")
					{
					    $o_counter++;
					}
					elsif($consensus_char eq "P" || $consensus_char eq "p")
					{
					    $p_counter++;
					}
					elsif($consensus_char eq "Q" || $consensus_char eq "q")
					{
					    $q_counter++;
					}
					elsif($consensus_char eq "R" || $consensus_char eq "r")
					{
					    $r_counter++;
					}
					elsif($consensus_char eq "S" || $consensus_char eq "s")
					{
					    $s_counter++;
					}
					elsif($consensus_char eq "T" || $consensus_char eq "t")
					{
					    $t_counter++;
					}
					elsif($consensus_char eq "U" || $consensus_char eq "u")
					{
					    $u_counter++;
					}
					elsif($consensus_char eq "V" || $consensus_char eq "v")
					{
					    $v_counter++;
					}
					elsif($consensus_char eq "W" || $consensus_char eq "w")
					{
					    $w_counter++;
					}
					elsif($consensus_char eq "X" || $consensus_char eq "x")
					{
					    $x_counter++;
					}
					elsif($consensus_char eq "Y" || $consensus_char eq "y")
					{
					    $y_counter++;
					}
					elsif($consensus_char eq "Z" || $consensus_char eq "z")
					{
					    $z_counter++;
					}
					elsif($consensus_char eq "-")
					{
					    $gap_counter++;
					}
				    }
				    
				    #...then each complete residue count is placed in its corresponding spot in a hash...
				    $Char_hash{$A} = $a_counter;
				    $Char_hash{$B} = $b_counter;
				    $Char_hash{$C} = $c_counter;
				    $Char_hash{$D} = $d_counter;
				    $Char_hash{$E} = $e_counter;
				    $Char_hash{$F} = $f_counter;
				    $Char_hash{$G} = $g_counter;
				    $Char_hash{$H} = $h_counter;
				    $Char_hash{$I} = $i_counter;
				    $Char_hash{$J} = $j_counter;
				    $Char_hash{$K} = $k_counter;
				    $Char_hash{$L} = $l_counter;
				    $Char_hash{$M} = $m_counter;
				    $Char_hash{$N} = $n_counter;
				    $Char_hash{$O} = $o_counter;
				    $Char_hash{$P} = $p_counter;
				    $Char_hash{$Q} = $q_counter;
				    $Char_hash{$R} = $r_counter;
				    $Char_hash{$S} = $s_counter;
				    $Char_hash{$T} = $t_counter;
				    $Char_hash{$U} = $u_counter;
				    $Char_hash{$V} = $v_counter;
				    $Char_hash{$W} = $w_counter;
				    $Char_hash{$X} = $x_counter;
				    $Char_hash{$Y} = $y_counter;
				    $Char_hash{$Z} = $z_counter;
				    $Char_hash{$gap} = $gap_counter;
				    

				    
				    
				    #...then, sort the hash keys by hash values (from largest to smallest) to determine the most frequently-occuring residue...
				    @sorted = sort{$Char_hash{$b} <=> $Char_hash{$a}} keys %Char_hash;
;
				    
				    if($Char_hash{$sorted[0]} == $Char_hash{$sorted[1]})
				    {#...which either produces a tie, so randomly select which to use...
					$random = rand;
					if($random < 0.5)
					{
					    $consensus_sequence .= $sorted[0];
					}
					else
					{
					    $consensus_sequence .= $sorted[1];
					}
				    }
				    else
				    {#...or puts the given residue at the beginning of the array, so simply add this element.	
				    $consensus_sequence .= $sorted[0]; #ultimately, this uses the most frequent residue at this position
				    
				    
				}
			}
			$gene_counter++;
			$gene_number = "protein_family_";
			$gene_number .= $gene_counter;
			print OUT ">$gene_number\n$consensus_sequence\n\n";
		}

		$consensus = '';
		$consensus_sequence = '';
		$gene_name = '';
		$protein_sequence = '';
		$consensus = '';
		undef %Protein_hash;
	}
	elsif(length($buf) < 4 && $buf !~ /\*/) #several blank lines exist between and within each alignment
	{
		$is_line_blank = "yes";
	}
	elsif($buf =~ /\w/ && $buf !~ /\(/) #indicates a line with protein name and partial sequence
	{		
	    $measure_space = $buf;
	    $offset = rindex($measure_space, " ");
	    ($gene_name, $temp_protein_sequence) = split(/\s+/, $buf);
	    $Protein_hash{$gene_name} .= $temp_protein_sequence;
	    $is_line_blank = "no";
		

	}
	elsif(length($buf) >= 4 || $buf =~ /\*/)
	{
		$raw_consensus = $buf;
		$temp_consensus = substr($raw_consensus, ($offset+1));
		$consensus .= $temp_consensus;
	}
	$iterator++;
	print "Currently on line $iterator of $in\n" if ($iterator % 100000 == 0);

}


$end_time = localtime;

print "\n\n$end_time\n";




