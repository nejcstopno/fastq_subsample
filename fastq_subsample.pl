#!/usr/bin/perl
use strict;
use warnings;

##################################################################################################
##################################################################################################
##												##
## Usage: perl fastq_subsample.pl fastq_file1 fastq_file2 genome_size coverage_desirable	##
##												##
##################################################################################################
##################################################################################################

# Public hash variables
my %hash_fasta1;
my %hash_quality1;
my %hash_fasta2;
my %hash_quality2;
my %hash_num_id;
my $total_fastq_nucleotides = 0;

# Save fastq file into 2 publics hash tables
sub fastq_to_hash()
{
	my $fastq_file = shift;
	open(FILE, $fastq_file);
	my @line;
	my $k = 1;
	my $seq_num = 1;
	$total_fastq_nucleotides = 0;
	while(<FILE>)
	{
		my $current_line = $_;
		chomp $current_line;
		if($k % 4 == 1) { $line[0] = $current_line; }
		elsif($k % 4 == 2) { $line[1] = $current_line; }
		elsif($k % 4 == 3) { $line[2] = $current_line; }
		else
		{
			$line[3] = $current_line;
			# forward fastq
			$hash_fasta1{$line[0]} = $line[1];
			$hash_quality1{$line[0]} = $line[3];
			# reverse fastq
			$hash_fasta2{$line[0]} = $line[1];
			$hash_quality2{$line[0]} = $line[3];
			# count nucleotides
			$total_fastq_nucleotides += length($line[1]);
			$hash_num_id{$seq_num} = $line[0];
			$seq_num++;
		}
		$k++;
	}
	close(FILE);
}

# print fastq hash
sub display_fastq_hash()
{
	foreach my $key (keys(%hash_fasta1))
	{
		print $key."\n".$hash_fasta1{$key}."\n".$hash_quality1{$key}."\n";
	}
}

# calculate coverage from genome size (in MB) hash_fastq
sub calculate_coverage()
{
        my $genome_size = shift;
        my $coverage;
        $coverage = int($total_fastq_nucleotides / ($genome_size * 1000000));
        return($coverage);
}

# calculate number of reads necessary to cover the desirable coverage
# get the average length of reads and return number of reads necessary to cover
# the desirable coverage
sub calculate_reads()
{
	my $genome_size = shift;
	my $total_reads = shift;
	my $nucleotides_to_complete_coverage = $genome_size * 1000000;
	my $necessary_reads = int($nucleotides_to_complete_coverage/($total_fastq_nucleotides/$total_reads));
	return($necessary_reads);
}

# generate an array with n no-duplicate random nums
sub rand_array()
{
        my $n = shift;
	my $total = shift;
        my @original_array = map { int } (1..$total);
        my @rand_array;
        my @sorted_array;
        my $size = scalar(@original_array);
        if( $n < $size && $n > 0 )
        {
                for (my $i = 0; $i <= $n; $i++)
                {
                        $size = $size -1;
                        my $index = int(rand($size));
                        $rand_array[$i] = $original_array[$index];
                        $original_array[$index] = $original_array[$size];
                }
                @sorted_array = sort {$a <=> $b} @rand_array;
                return(@sorted_array);
        }
        else
        {
                return(@original_array);
        }
}

# save sub hashes in 2 new files
sub sub_fastq_hash()
{
	my $num_of_reads = shift;
	my $total_reads = shift;
	my $file_forward_fastq = shift;
	my $file_reverse_fastq = shift;
	my @array_random_index = &rand_array($num_of_reads, $total_reads);
	open(FILE1, ">", $file_forward_fastq) or die $!;
	open(FILE2, ">", $file_reverse_fastq) or die $!;
	foreach (@array_random_index)
	{
		my $key = $hash_num_id{$_};
		print FILE1 $key."\n".$hash_fasta1{$key}."\n+\n".$hash_quality1{$key}."\n";
		print FILE2 $key."\n".$hash_fasta2{$key}."\n+\n".$hash_quality2{$key}."\n";
	}
	close(FILE1);
	close(FILE2);
}

# main
my $fastq_file1 = $ARGV[0];
my $fastq_file2 = $ARGV[1];
my $genome_size = $ARGV[2];
my $desirable_coverage = $ARGV[3];

&fastq_to_hash($fastq_file1);

my $coverage = &calculate_coverage($genome_size);
my $num_of_reads = &calculate_reads($genome_size,$desirable_coverage);
my $total_reads = scalar keys %hash_fasta1;
#my @array_rand = &calculate_reads($num_of_reads, $total_reads);
print "Estimated coverage: ".$coverage."\n";
print "Desirable coverage: ".$desirable_coverage."\n";
my @array = &rand_array($num_of_reads, $total_reads);
print join(", ", @array);
#print &rand_array($num_of_reads, $total_reads);
&sub_fastq_hash($num_of_reads, $total_reads, "sub_".$fastq_file1, "sub_".$fastq_file2);
print "See files sub_".$fastq_file1."sub_".$fastq_file2."\n";
