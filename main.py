from flask import Flask, render_template, send_file, request

app = Flask(__name__)


@app.route('/')
def main():
    return render_template('index.html')


@app.route('/database/chipseq_index')
def chipseq():
    return send_file('sample_data/ChIPseq.index')

@app.route('/database/chipseq_sample')
def chipseq_sample():
    return send_file(F'sample_data/ChIP-seq/{request.args.get("sample")}')

@app.route('/database/null_model')
def null_model():
    return send_file('sample_data/Null/CommonSNP_MAFmatch')

# TODO: Will send dbSNP data file, but currently too big to fit in GitHub.