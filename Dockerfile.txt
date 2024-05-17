# Use an official Python runtime as a parent image
FROM continuumio/miniconda3

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# Install system dependencies
RUN apt-get update && \
    apt-get install -y libxext6

# Set LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib

# Create and activate Conda environment
RUN conda create -y --name for_serotoninAI
RUN echo "source activate for_serotoninAI" > ~/.bashrc
ENV PATH /opt/conda/envs/for_serotoninAI/bin:$PATH

# Install RDKit and other dependencies
RUN conda install -y -c conda-forge rdkit==2023.03.2
RUN conda install -y -c anaconda scikit-learn==1.2.2
RUN pip install ipython
RUN pip install wordcloud==1.9.1.1
RUN conda install -y -c conda-forge mljar-supervised
RUN conda install -y -c anaconda numpy==1.26.4 
RUN conda install -y -c anaconda pandas==1.5.3
RUN conda install -y -c anaconda seaborn==0.12.2
RUN conda install -y -c conda-forge matplotlib==3.4.3
RUN conda install -y -c conda-forge mordred==1.2.0
RUN pip install scipy==1.11.4
RUN pip install streamlit==1.27.0
RUN pip install streamlit-option-menu==0.3.6
RUN pip install streamlit-ketcher==0.0.1
RUN pip install supervised

# Set the working directory in the container
WORKDIR /app

# Copy the Streamlit app file into the container at /app
COPY app_streamlit_SerotoninAI.py /app/

# Expose the port the app runs on
EXPOSE 8501

# Set entrypoint to run Streamlit app
ENTRYPOINT ["streamlit", "run", "app_streamlit_SerotoninAI.py"]
