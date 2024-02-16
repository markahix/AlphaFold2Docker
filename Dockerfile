FROM continuumio/miniconda3:latest
# Base Flask Application Setup
RUN pip install flask
RUN apt-get update
RUN apt-get install -y --no-install-recommends libatlas-base-dev gfortran nginx supervisor
RUN pip install  --no-binary pyuwsgi pyuwsgi

# Application-Specific Setup
RUN pip install --no-warn-conflicts 'colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold'
RUN pip install --upgrade dm-haiku
RUN ln -s /opt/conda/lib/python3.11/site-packages/colabfold colabfold
RUN ln -s /opt/conda/lib/python3.11/site-packages/alphafold alphafold
RUN pip install --upgrade "jax[cuda11_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
RUN pip install ml_dtypes==0.2.0
RUN sed -i 's/weights = jax.nn.softmax(logits)/logits=jnp.clip(logits,-1e8,1e8);weights=jax.nn.softmax(logits)/g' alphafold/model/modules.py
RUN apt install -y zip

COPY Application/ /app
WORKDIR /app

EXPOSE 8543
# CMD ["uwsgi","--http","0.0.0.0:8543","--master","-p","4","-w","app:app"]
CMD ["flask","run","--host=0.0.0.0","--port=8543"]

