from flask import Flask, render_template, request
from werkzeug.utils import secure_filename
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

app = Flask(__name__)

app.config["UPLOAD_FOLDER"] = "data/"
E_VALUE_THRESH = 0.04

@app.route('/')
def index():
	return render_template("index.html")

@app.route('/results', methods=["POST"])
def form():
	#takes text sequence and puts it in to input.fasta
	if(request.form["submitButton"] == "text"):
		filename = app.config['UPLOAD_FOLDER'] + "input.fasta"
		sequence = request.form.get("sequence")
		fasta_file = open(filename,"a+")
		fasta_file.truncate(0)
		fasta_file.write(sequence)
		fasta_file.close()
	#check file input FASTA format
	elif(request.form["submitButton"] == "file"):
		fileInput = request.files["file"]
		filename = app.config['UPLOAD_FOLDER'] + secure_filename(fileInput.filename)
		if('.' in filename and filename.rsplit('.', 1)[1] != "fasta"):
			error_statement = "File not in FASTA format"
			return render_template("index.html", error_statement=error_statement)
		fileInput.save(filename)
	
	#check if input is in FASTA format	
	if(is_fasta(filename) == False):
			error_statement = "Data not in FASTA format"
			return render_template("index.html", error_statement=error_statement)

	#parse result and read alignments into array to display on table
	record = SeqIO.read(filename, format="fasta")
	result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
	blast_record = NCBIXML.read(result_handle)
	alignArray = []
	for alignment in blast_record.alignments:
	    for hsp in alignment.hsps:
	        if hsp.expect < E_VALUE_THRESH:
	            align = AlignObject(alignment.title, alignment.length, hsp.score, hsp.gaps, hsp.expect)
	            alignArray.append(align)

	return render_template("results.html", alignArray=alignArray)

#Checks if input file is in FASTA format
def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

#object to store alignment data
class AlignObject:
	def __init__(obj, title, length, score, gaps, expect):
		obj.title = title
		obj.length = length
		obj.score = score
		obj.gaps = gaps
		obj.expect = expect