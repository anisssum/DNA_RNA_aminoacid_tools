import os
import unitest
import pytest
from DNA_RNA_aminoacid_tools import SHORT_CODE, run_genescan, filter_fastq
from bio_files_processor import FastaRecord

def test_sort_code_len():
	target = 44
	result = len(SHORT_CODE)
	assert target == result

def test_short_code_content():
	target = [
		'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U', 'O',
    		'a', 'r', 'n', 'd', 'c', 'e', 'q', 'g', 'h', 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v', 'u', 'o'
	]
	result = SHORT_CODE
	assert target == result

def test_run_genescan_error():
	with pytest.raises(ValueError("Neither sequence nor sequence_file provided")):
		run_genescan()

@pytest.fixture
def tmp_file():
    file_path = './test/test_filtered.fastq'
    yield file_path
    if os.path.exists(file_path):
        os.remove(file_path)

def test_filter_fastq_exists(tmp_file):
    filter_fastq('./test_data/test.fastq', output_filename = tmp_file)
    assert os.path.exists(tmp_file)

def test_show_length_correct():
	inp = 'qchimfwr'
	target = 8
	result = show_length(inp)
	assert target == result

def test_FastaRecord_repr():
    fasta_record = FastaRecord(id = 'test_id', description = 'test_description', sequence = 'ATCGATCGAT')
    expected_repr = "FastaRecord(id='test_id', description='test_description', sequence='ATCGATCGAT...')"
    assert repr(fasta_record) == expected_repr

def test_molecular_weight_correct():
	inp = 'omg'
	target = 443.55
	result = molecular_weight(inp)
	assert target == result