use LWP::Simple;
$query = 'eukaryota%5Bporgn%5D+NOT+(metazoa%5Bporgn%5D+OR+fungi%5Bporgn%5D+OR+viridiplantae%5Bporgn%5D)+AND+srcdb_refseq%5Bprop%5D+AND+biomol_genomic%5Bprop%5D+AND+(NC_000000%3ANC_999999%5Bpacc%5D+OR+AC_000000%3AAC_999999%5Bpacc%5D+OR+(NT_000001%3ANT_999999999%5Bpacc%5D+AND+(%22chromosome+2L%22+OR+%22chromosome+2R%22+OR+%22chromosome+3L%22+OR+%22chromosome+3R%22)))';
#assemble the esearch URL
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "esearch.fcgi?db=nuccore&term=${query}&usehistory=y";

#post the esearch URL
$output = get($url);

#parse WebEnv, QueryKey and Count (# records retrieved)
$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
$count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

#open output file for writing
my $filename = 'other_eukaryota_chrom.fa';
open(my $fh, '>', $filename) or die "Can't open file!\n";

#retrieve data in batches of 500
$retmax = 500;
for ($retstart = 0; $retstart < $count; $retstart += $retmax) {
        $efetch_url = $base ."efetch.fcgi?db=nucleotide&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
        $efetch_out = get($efetch_url);
        print $fh "$efetch_out";
}
close $fh;
