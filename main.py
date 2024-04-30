from flask import Flask, render_template, send_file, request
#import processing_file_name
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

# @app.route('/database/null_model')
# def null_model():
#     return send_file('sample_data/Null/CommonSNP_MAFmatch')

@app.route('/upload', methods=['POST'])
def upload_file():
    uploaded_file = request.files['file']
    print("Uploaded file:", uploaded_file)
    if uploaded_file.filename != '':
        uploaded_file.save("test")
    #results = processing_file_name.some_function(uploaded_file)
    #return results
    

# TODO: Will send dbSNP data file, but currently too big to fit in GitHub.