# SerotoninAI_streamlit

The web page application for affinity prediction for serotonergic targets and properties (Blood-brain barrier penetration and Human intestinal absorption)

![Logo SerotoninAI](https://github.com/nczub/SerotoninAI_streamlit/blob/main/serotoninAI_homepage.png)

## Online Access
You can find the application at the following link:
[https://serotoninai.streamlit.app/](https://serotoninai.streamlit.app/)

## Docker with SerotoninAI
https://hub.docker.com/r/serotoninai/serotoninai-app

## Installation for Local Development

To run the application locally, follow these steps:
1. **Clone the repository**
2. **Install dependencies**:

The needed packages are in file enviroment.txt

During installation you create conda environment 'for_serotoninAI'

3. **Activate environment**
   
In the console activate conda environment:

```bash
$ conda activate for_serotoninAI
```

4. Now, you can run the application:

```bash   
$ streamlit run app_streamlit_SerotoninAI.py
```
App should open in the browser or it will be available at 'http://localhost:8501'.

5. Finally, have fun and test my app!


## Batch mode

In batch mode, you can calculate predictions for multiple molecules. The online version of SerotoninAI has limitations based on Streamlit Cloud. The local app is much better to use for larger database.

Please, remember to upload CSV file with the column names 'smiles', because based on this system will predict affinity or property.

## Author

I'm Natalia Łapińska (maiden name Czub) and I'm the author of SerotoninAI

- GitHub: https://github.com/nczub

- LinkedIn: https://www.linkedin.com/in/natalia-czub/

## Articles
- SerotoninAI: Serotonergic System Focused, Artificial Intelligence-Based Application for Drug Discovery

*Natalia Łapińska, Adam Pacławski, Jakub Szlęk, and Aleksander Mendyk*
https://doi.org/10.1021/acs.jcim.3c01517

- Integrated QSAR Models for Prediction of Serotonergic Activity: Machine Learning Unveiling Activity and Selectivity Patterns of Molecular Descriptors

*Natalia Łapińska, Adam Pacławski, Jakub Szlęk, and Aleksander Mendyk*
https://doi.org/10.3390/pharmaceutics16030349

## License

This project is available under the GNU General Public License v3.0 (GPL-3.0).
