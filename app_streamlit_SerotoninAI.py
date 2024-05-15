#Natalia Łapińska

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.


#load packages
import streamlit as st
import requests
from streamlit_option_menu import option_menu
from rdkit import Chem
from rdkit.Chem import Draw
from mordred import Calculator, descriptors
from supervised.automl import AutoML
import pandas as pd
import numpy as np
import time
import hashlib
import json
import html
from streamlit.components.v1 import html
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from PIL import Image
from streamlit_ketcher import st_ketcher
import matplotlib.pyplot as plt
import smtplib
from email.message import EmailMessage

st.set_page_config(page_title="SerotoninAI")

#background of webpage
page_bg_img = f"""
<style>
[data-testid="stAppViewContainer"] > .main {{
background-image: url("https://raw.githubusercontent.com/nczub/APP/2475af8d4e836424e301f5829c83cbbf8632d2dc/background_14_07_2023.svg");
background-size: cover;
background-position: top;
#background-repeat: repeat;
background-attachment: local;
}}
[data-testid="stHeader"] {{
background: rgba(0,0,0,0);

}}
</style>
"""
st.markdown(page_bg_img, unsafe_allow_html=True)



#loading models
#5-HT1A
best_model_path_5HT1A = 'FINAL_QSAR_MODELS/model_5-HT1A/mljar_AutoML_Compete_2023_07_24_16_07_19' 
best_model_5HT1A = AutoML(best_model_path_5HT1A)

#5-HT1B
best_model_path_5HT1B = 'FINAL_QSAR_MODELS/model_5-HT1B/mljar_AutoML_Compete_2023_07_20_23_26_40' 
best_model_5HT1B = AutoML(best_model_path_5HT1B)

#5-HT1D
best_model_path_5HT1D = 'FINAL_QSAR_MODELS/model_5-HT1D/mljar_AutoML_Compete_2023_07_22_05_05_38' 
best_model_5HT1D = AutoML(best_model_path_5HT1D)

#5-HT2A
best_model_path_5HT2A = 'FINAL_QSAR_MODELS/model_5-HT2A/mljar_AutoML_Compete_2023_07_24_15_24_54' 
best_model_5HT2A = AutoML(best_model_path_5HT2A)

#5-HT2B
best_model_path_5HT2B = 'FINAL_QSAR_MODELS/model_5-HT2B/mljar_AutoML_Compete_2023_07_21_12_27_51' 
best_model_5HT2B = AutoML(best_model_path_5HT2B)

#5-HT2C
best_model_path_5HT2C = 'FINAL_QSAR_MODELS/model_5-HT2C/mljar_AutoML_Compete_2023_07_24_15_35_02' 
best_model_5HT2C = AutoML(best_model_path_5HT2C)

#5-HT3
best_model_path_5HT3 = 'FINAL_QSAR_MODELS/model_5-HT3/mljar_AutoML_Compete_2023_07_21_16_02_20' 
best_model_5HT3 = AutoML(best_model_path_5HT3)

#5-HT4
best_model_path_5HT4 = 'FINAL_QSAR_MODELS/model_5-HT4/mljar_AutoML_Compete_2023_11_08_16_50_53'
best_model_5HT4 = AutoML(best_model_path_5HT4)

#5-HT5A
best_model_path_5HT5A = 'FINAL_QSAR_MODELS/model_5-HT5A/mljar_AutoML_Compete_2023_07_21_17_43_12' 
best_model_5HT5A = AutoML(best_model_path_5HT5A)

#5-HT6
best_model_path_5HT6 = 'FINAL_QSAR_MODELS/model_5-HT6/mljar_AutoML_Compete_2023_07_24_14_35_11' 
best_model_5HT6 = AutoML(best_model_path_5HT6)

#5-HT7
best_model_path_5HT7 = 'FINAL_QSAR_MODELS/model_5-HT7/mljar_AutoML_Compete_2023_07_25_10_01_44' 
best_model_5HT7 = AutoML(best_model_path_5HT7)

#SERT
best_model_path_sert = 'FINAL_QSAR_MODELS/model_SERT/mljar_AutoML_Compete_2023_07_25_10_22_08' 
best_model_sert = AutoML(best_model_path_sert)

#HIA models
#regression
best_model_path_hia_regression = 'FINAL_QSPR_MODELS/regression_HIA/mljar_AutoML_Compete_2022_11_15_13_28_31' 
best_model_hia_regression = AutoML(best_model_path_hia_regression)

#classification
best_model_path_hia_classification = 'FINAL_QSPR_MODELS/classification_HIA/mljar_AutoML_Compete_2022_09_26_13_38_07' 
best_model_hia_classification = AutoML(best_model_path_hia_classification)

#BBB model
best_model_path_BBB = 'FINAL_QSPR_MODELS/BBB_model/mljar_AutoML_Compete_2023_08_17_14_17_34'
best_model_BBB = AutoML(best_model_path_BBB)

#binary model for active and inactive compounds - serotonergic activity
binary_model_path = 'FINAL_QSAR_MODELS/serotonergic_activity'

#selectivity path
selective_model_path = 'FINAL_QSAR_MODELS/selectivity'



#for mordred descriptors calculations
descriptors_for_QSAR = pd.read_table("descriptors_QSAR.txt", header = None)[0][0]
calc = Calculator(descriptors, ignore_3D=True)
calc.descriptors = [d for d in calc.descriptors if str(d) in descriptors_for_QSAR] 

#clearning SMILES input button
def clear_text():
    st.session_state["text"] = ""

custom_css = """
<style>
    :root {
        font-size: 20px;
        text-align: justify;
    }
    .text-second-title {
        font-size: 40px;
        text-align: left;
        color: #525354;
    }

    @keyframes text-gradient-title {
        0% { color: grey; }
        50% { color: #91b3bd; }
        100% { color: grey; }
    }

    .text-gradient-title {
        position: sticky;
        top: 0px;
        animation: text-gradient-title 4s ease-in-out infinite;
        font-size: 130px;
        text-align: center;
        text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.8);
        font-style: italic;
    }
    .centered-image {
        display: flex;
        justify-content: center;
    }
    </style>

</style>
"""

st.markdown(custom_css, unsafe_allow_html=True)
#main title
st.markdown('<h1 class="text-gradient-title">SerotoninAI</h1>', unsafe_allow_html=True)


success_style = """
    background-color: #b0d1e0;
    color: #525354;
    border-radius: 10px;
    padding: 10px;
    width: 80px;
    fontSize: 25px;
    animation-name: fadeOut;
    animation-duration: 5s;
"""

#footer
footer_style = """
    position: fixed;
    left: 0;
    z-index: 3;
    bottom: 0;
    width: 100%;
    color: #525354;
    font-style: italic;
    text-align: left;
    padding: 10px;
    font-size: 16px;
"""
st.markdown(
    f'<div style="{footer_style}">Copyright (C) 2024 Natalia Łapińska (Czub)</div>',
    unsafe_allow_html=True
)

if st.session_state.get('switch_button', False):
    st.session_state['menu_option'] = (st.session_state.get('menu_option', 0) + 1) % 2
    manual_select = st.session_state['menu_option']
else:
    manual_select = None

selected = option_menu(None, ["Home", "5-HT receptors", "SERT", "Batch calculation", "HIA", "BBB", "Serotonergic activity", "Selectivity", "Q&A",  "Contact"],
                        icons=['house', "layers", "layers-fill",'cloud-upload', "capsule-pill",  "capsule",  "activity", "fingerprint", "question-diamond",'envelope-at'],
                        orientation="horizontal", manual_select=manual_select, key='menu_20', default_index = 0,
                        styles={
        "container": {"padding": "21!important", "background-color":"#b4bbbf", "width": "auto"},
        "icon": {"color": "#4e5152", "font-size": "21px", "text-align" : "center"}, 
        "nav-link": {"font-size": "25px", "text-align": "center", "margin":"5px", "--hover-color": "#757473", "font-color":"#0a0a0a"},
        "nav-link-selected": {"background-color": "#5d93a3"},
        })
if selected == "Home":
    st.markdown("<h1 style='text-align: center; fontSize: 30px; font-style: italic; color: grey;'>Serotonin Science Revolution: Decipher Molecule's Affinity and Properties for Future Breakthrough Therapy!</h1>", unsafe_allow_html=True)
    st.write('')
    st.markdown("<h1 style='text-align: center; fontSize: 23px; color: grey; font-weight: normal;'>Discover instant affinity predictions with our cutting-edge app, propelling your breakthroughs to new heights. Join the revolution shaping serotonergic studies and leave your mark on scientific progress!</h1>", unsafe_allow_html=True)
    st.write("---")
    st.subheader("Instruction")
    st.write("This application provides the ability to predict the ligand affinity for a specific serotonergic target, as well as human intestinal absorption, which is important for orally administered drugs. Moreover, the app gives you predictions for blood-brain barrier penetration.")
    st.write("In the menu above you can find the main sections:")
    st.write(":small_blue_diamond: **5-HT receptors**")
    st.write("The possibility of calculating the affinity for specific serotonergic receptors: 5-HT1A, 5-HT1B, 5-HT1D, 5-HT2A, 5-HT2B, 5-HT2C, 5-HT3, 5-HT4, 5-HT5A, 5-HT6 and 5-HT7.")
    st.write("As you probably have noticed this app does not calculate predictions for all serotonin receptors. The reason behind this is limited information about ligands and their affinity to 5-HT1E, 5-HT1F, and 5-HT5B receptors. We decided not to provide models which were built on a dataset smaller that 500 molecules. If Quantitative structure–activity relationship models (QSAR models) are based on a minor database, such a model could be overfitted to  the training set. As these databases grow, the application will be updated.")
    st.write('')
    st.write(":small_blue_diamond: **SERT**")
    st.write("Among serotonergic receptors, the affinity for the serotonin transporter can be predicted on this subpage.")
    st.write('')
    st.write(":small_blue_diamond: **Batch calculation**")
    st.write("Sometimes you need to work more efficiently, so we created a batch mode. It gives you the ability to predict affinity for selected targets and for more than one molecule at a time.")
    st.write('')
    st.write(":small_blue_diamond: **HIA**")
    st.write("Human intestinal absorption (HIA) is a property of drugs based on how efficiently they pass through the intestinal wall. Drugs with high HIA can be administered orally.")
    st.write("On this subpage you can calculate the prediction of a specific class - low or high absorption for a single or multiple molecules.")
    st.write("Due to the fact that we used only molecules with serotonergic activity to train the artificial intelligence-based system, we do not recommend using it with other biological ligand groups. This may result in incorrect predictions.")
    st.write('')
    st.write(':small_blue_diamond: **BBB**')
    st.write('The blood-brain barrier is a barrier between blood vessels and nervous tissue, designed to protect the nervous system against harmful factors, and to enable the selective transport of substances from the blood to the cerebrospinal fluid.')
    st.write('On this page you can get BBB penetration prediction. This is the last element to predict the key actions and properties of serotonergic compounds.')
    st.write('')
    st.write(":small_blue_diamond: **Q&A**")
    st.write("The question and answer subpage allows you to review possible problems when using this application.")
    st.write('')
    st.write(":small_blue_diamond: **Contact**")
    st.write("The last subpage contains information on how to contact the author of **SerotoninAI**.")
    st.write('')
    st.write('---')
    st.subheader("References")
    st.write("The provided information and graphics are made possible thanks to the following sources:")
    st.write(":small_blue_diamond: Sharp T, Barnes NM. Central 5-HT receptors and their function; present and future. Neuropharmacology. 2020;177:108155. doi:10.1016/j.neuropharm.2020.108155.")
    st.write(":small_blue_diamond: Rudnick G, Sandtner W. Serotonin transport in the 21st century. J Gen Physiol. 2019;151(11):1248-1264. doi:10.1085/jgp.201812066.")
    st.write(":small_blue_diamond: Newman-Tancredi A, Depoortère RY, Kleven MS, Kołaczkowski M, Zimmer L. Translating biased agonists from molecules to medications: Serotonin 5-HT1A receptor functional selectivity for CNS disorders. Pharmacol Ther. 2022;229:107937. doi:10.1016/j.pharmthera.2021.107937.")
    st.write(":small_blue_diamond: Czub N, Szlęk J, Pacławski A, Klimończyk K, Puccetti M, Mendyk A. Artificial Intelligence-Based Quantitative Structure-Property Relationship Model for Predicting Human Intestinal Absorption of Compounds with Serotonergic Activity. Mol Pharm. 2023;20(5):2545-2555. doi:10.1021/acs.molpharmaceut.2c01117.")
    st.write(":small_blue_diamond: Sivandzade F, Cucullo L. In-vitro blood-brain barrier modeling: A review of modern and fast-advancing technologies. J Cereb Blood Flow Metab. 2018;38(10):1667-1681. doi:10.1177/0271678X18788769.")
    st.write(":small_blue_diamond: Plonska, A.; Plonski, P. MLJAR: State-of-the-Art Automated Machine Learning Framework for Tabular Data, version 0.10.3. https://github.com/mljar/mljar-supervised.")
    st.write(":small_blue_diamond: The graphics were created using Canva, a user-friendly design tool. https://www.canva.com/")
    st.write('---')
    st.subheader('Acknowledgments')
    st.write('💎 I would like to thank my team **Aleksander Mendyk**, **Jakub Szlęk** and **Adam Pacławski** for their invaluable contribution. Your collective knowledge laid the foundation for my research journey. Your unwavering support and dedication have been instrumental. Your collaborative spirit in creating models for 5-HT1A and HIA has been remarkable. I also appreciate your insightful suggestions that have enriched the refinement and evolution of ***SerotoninAI***.')
    st.write('💎 In addition, I want to thank **Klaudia Klimończyk** as a co-author of the HIA module.')
    st.write('💎 My thanks extend to **Justyna Srebro** for her support and graphic suggestions.')
    st.write('💎 Finally, I would like to thank **Paweł Łapiński** for his encouragement and language correction.')
    st.write('---')
    st.subheader('License')
    st.write('GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007')
    st.write('Copyright (C) 2024 Natalia Łapińska (Czub)')
    st.write('This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. https://www.gnu.org/licenses/gpl-3.0.html')
    st.write('DISCLAIMER OF LIABILITY')
    st.write("The author of this software shall not be liable for any special, incidental, consequential, or indirect damages resulting from the use, misuse, or inability to use this software, including but not limited to, damages for loss of profits, business interruption, or loss of data. The software is provided 'as is' and the author make no warranties, either express or implied, regarding the software's fitness for a particular purpose or its accuracy and reliability.")


elif selected == "5-HT receptors":
    calc = Calculator(descriptors, ignore_3D=True)
    calc.descriptors = [d for d in calc.descriptors if str(d) in descriptors_for_QSAR] 
    st.markdown('<h1 class="text-second-title">Choose specific receptor</h1>', unsafe_allow_html=True)
    sub_option = st.selectbox("", ["5-HT1A receptor", "5-HT1B receptor","5-HT1D receptor",
                                                  "5-HT2A receptor","5-HT2B receptor","5-HT2C receptor",
                                                  "5-HT3 receptor", "5-HT4 receptor", "5-HT5A receptor","5-HT6 receptor","5-HT7 receptor"])
    if sub_option == "5-HT1A receptor":
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            predictions = best_model_5HT1A.predict(descriptors_value_df)
                        st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                        prediction_float = round(float(predictions), 3)
                        st.write("pKi value for 5-HT1A serotonin receptor: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['BCUTc-1l', 'MAXdO', 'MAXaaaC', 'PEOE_VSA9', 'SMR_VSA3', 'SssCH2', 
                                 'AATS4i', 'SpAbs_DzZ', 'AATS6dv', 'VSA_EState5']
                        min_values = {'BCUTc-1l': -0.7373620145398293, 'MAXdO': 9.099537037037038, 'MAXaaaC': -0.1340347251287001,
                        'PEOE_VSA9': 0.0, 'SMR_VSA3': 0.0, 'SssCH2': -0.4661989795918362, 'AATS4i': 147.56632501478904,
                        'SpAbs_DzZ': 42.05895519992669, 'AATS6dv': 0.1538461538461538, 'VSA_EState5': -7.181078813682171}
                        max_values = {'BCUTc-1l': -0.292146392415571, 'MAXdO': 14.629783178882205, 'MAXaaaC': 1.5381250000000002,
                        'PEOE_VSA9': 78.6625871223764, 'SMR_VSA3': 39.8909626482546, 'SssCH2': 23.225582666887828,
                        'AATS4i': 175.1107976137481, 'SpAbs_DzZ': 1265.278990098867, 'AATS6dv': 7.298507462686567,
                        'VSA_EState5': 8.521368790538302}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)
        agree_draw_smiles = st.checkbox("Draw chemical structure")
        if agree_draw_smiles:
            smile_code = st_ketcher()
            st.markdown(f"SMILES: {smile_code}")
        st.write("---")
        st.subheader("5-HT1A in pharmacology")
        st.write("5-HT1A receptor is the best known serotonergic receptor and is abundantly present in cortical and limbic areas.")
        st.write("Agonists of the 5-HT1A receptor induce a variety of behavioral and physiological effects, including altered motor function, body temperature, and neuroendocrine activity. Furthermore, there is a wealth of preclinical and clinical evidence suggesting that 5-HT1A receptor agonists possess antidepressant and anxiolytic properties, and they impact various cognitive domains relevant to schizophrenia symptoms. The numerous and diverse behavioral effects of 5-HT1A receptor agonists likely reflect the actions of these receptors in multiple forebrain and midbrain regions. For instance, the activation of 5-HT1A autoreceptors has been associated with anxiolytic effects, while the activation of postsynaptic 5-HT1A receptors is linked to antidepressant effects.")
        col1, col2 = st.columns(2)
        with col1:
            st.write("In the 1990s the 5-HT1A agonist buspirone was approved for clinical use as an anxiolytic, and structurally similar compounds such as gepirone and tandospirone demonstrated efficacy as antidepressants in clinical trials. Novel antidepressant medications, like vilazodone and vortioxetine, which combine high affinity for the 5-HT transporter and the 5-HT1A receptor (as well as several other 5-HT receptors in the case of vortioxetine), have been developed.")
        with col2:
            image = Image.open('images_app/5-HT1A.png')
            st.image(image, width=360)
        st.write("Within the 5-HT1A receptor ligands, we can distinguish so-called biased-agonists, which are characterized by their ability to selectively activate specific signaling pathways. It is this unique feature of biased-agonists that allows us to differentiate their effects on the cell compared to traditional agonists.")
        st.write("NLX-112 (formerly known as F13640) is considered a potential treatment for Parkinson's disease patients suffering from L-DOPA-induced dyskinesia. NLX-112 is a full agonist but shows preferential activation of ERK signaling and activates 5-HT1A autoreceptors to inhibit the firing of 5-HT cells.")
        st.write("Another biased agonist of the 5-HT1A receptor, NLX-101 (F15599), has been reported to preferentially activate postsynaptic 5-HT1A receptors and exhibits antidepressant, pro-cognitive, and neuroprotective properties in animal models.")        
    elif sub_option == "5-HT1B receptor":
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            predictions = best_model_5HT1B.predict(descriptors_value_df)
                        st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                        prediction_float = round(float(predictions), 3)
                        st.write("pKi value for 5-HT1B serotonin receptor: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['MINaaaC', 'SMR_VSA6', 'SlogP_VSA1', 'NsssN', 'VSA_EState4', 
                                 'SaaaC', 'ABC', 'BCUTZ-1l', 'SlogP_VSA2', 'BCUTse-1l']
                        min_values = {'MINaaaC': 0.0430998207874111, 'SMR_VSA6': 0.0, 'SlogP_VSA1': 0.0, 'NsssN': 0, 'VSA_EState4': -3.877327289138031, 'SaaaC': 0.0,
                        'ABC': 8.485281374238573, 'BCUTZ-1l': 0.9918795477066924, 'SlogP_VSA2': 4.9839785209472085, 'BCUTse-1l': 2.3910240890446155}
                        max_values = {'MINaaaC': 1.4785704837490554, 'SMR_VSA6': 115.64600422703208, 'SlogP_VSA1': 31.57463806993713, 'NsssN': 6,
                        'VSA_EState4': 18.322880995240126, 'SaaaC': 4.997746037925896, 'ABC': 48.90440618467739, 'BCUTZ-1l': 5.772390550568072, 'SlogP_VSA2': 137.2601831474361,
                        'BCUTse-1l': 2.498598365922732}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")                             
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)
        agree_draw_smiles = st.checkbox("Draw chemical structure")
        if agree_draw_smiles:
            smile_code = st_ketcher()
            st.markdown(f"SMILES: {smile_code}")                  
        st.write("---")
        col1, col2 = st.columns([2,1])
        with col1:
            st.subheader("5-HT1B in pharmacology")
            st.write("The 5-HT1B receptors play a prominent role in the basal ganglia, mesolimbic pathways, and the trigeminal nucleus. However, their significance in neuropharmacology should not be underestimated. Currently, 5-HT1B/D agonists known as triptans are considered the gold standard for acute migraine treatment.")
        with col2:
            image = Image.open('images_app/5-HT1B.png')
            st.image(image, width=250)
        st.write("These receptors hold substantial relevance beyond migraine treatment. They are actively studied in the context of depression, anxiety, aggression, and impulse control. Both preclinical and clinical investigations are shedding light on their functions. Notably, conditional tissue-specific and time-dependent 5-HT1B receptor knockout studies in mice contribute to understanding the roles of presynaptic and postsynaptic 5-HT1B receptors in these behaviors.")
        st.write("Despite the availability of selective 5-HT1B receptor antagonists such as SB224289, progress in 5-HT1B receptor research has been hampered by the lack of highly selective and brain-penetrating 5-HT1B agonists. Nevertheless, the study of 5-HT1B receptors remains crucial in the field of neuropharmacology, particularly because of their involvement in migraine treatment and a broader understanding of various neuropsychiatric conditions.")

    elif sub_option == "5-HT1D receptor":
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            predictions = best_model_5HT1D.predict(descriptors_value_df)
                        st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                        prediction_float = round(float(predictions), 3)
                        st.write("pKi value for 5-HT1D serotonin receptor: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['BCUTse-1l', 'BCUTZ-1l', 'BalabanJ', 'MAXaaCH', 'n10FARing', 'MAXaasC', 'SlogP_VSA2', 'BCUTi-1l', 'JGI3', 'EState_VSA8']
                        min_values = {'BCUTse-1l': 2.39096093838768, 'BCUTZ-1l': 0.9918795477066924, 'BalabanJ': 4.133332842407474e-07, 'MAXaaCH': 1.0679894687443388,
                        'n10FARing': 0, 'MAXaasC': 0.1587137776339895, 'SlogP_VSA2': 4.9839785209472085,
                        'BCUTi-1l': 10.300242888872871, 'JGI3': 0.0166666666666666, 'EState_VSA8': 0.0}
                        max_values = {'BCUTse-1l': 2.498598365922729, 'BCUTZ-1l': 5.772390550568072, 'BalabanJ': 2.622001753479514, 'MAXaaCH': 2.510647543251294,
                        'n10FARing': 2, 'MAXaasC': 1.6150434618291762, 'SlogP_VSA2': 111.60791765215704, 'BCUTi-1l': 10.984220997218516,
                        'JGI3': 0.0852272727272727, 'EState_VSA8': 108.07472392106264}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)
        agree_draw_smiles = st.checkbox("Draw chemical structure")
        if agree_draw_smiles:
            smile_code = st_ketcher()
            st.markdown(f"SMILES: {smile_code}")
        st.write("---")
        st.subheader("5-HT1D in pharmacology")
        col1, col2 = st.columns([1,3])
        with col1:
            image = Image.open('images_app/5-HT1D.png')
            st.image(image, width=200)
        with col2:
            st.write("The 5-HT1D receptors are abundantly present in the trigeminal nucleus. However, their significance in neuropharmacology should not be underestimated. Currently, triptans, which are agonists of both the 5-HT1B and 5-HT1D receptors, are considered the gold standard for acute migraine treatment. Nevertheless, the relative importance of 5-HT1B compared to 5-HT1D in this context remains uncertain.")
            st.write("The exploration of 5-HT1D receptors holds promise for the development of novel treatment approaches in migraine and other related disorders.")
        

    elif sub_option == "5-HT2A receptor":
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            predictions = best_model_5HT2A.predict(descriptors_value_df)
                        st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                        prediction_float = round(float(predictions), 3)
                        st.write("pKi value for 5-HT2A serotonin receptor: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ["nFRing", "VR1_A", "n6ARing", "SMR_VSA7", "nBase", "BCUTs-1h", "BCUTdv-1l", "BCUTZ-1l", "GATS5i", "BCUTdv-1h"]
                        min_values = {'nFRing': 0, 'VR1_A': 30.671729502867567, 'n6ARing': 0, 'SMR_VSA7': 0.0, 'nBase': 0, 'BCUTs-1h': 2.164333480214669,
                        'BCUTdv-1l': 0.1514613373305368, 'BCUTZ-1l': 5.636495147426821, 'GATS5i': 0.375280485256885, 'BCUTdv-1h': 4.048936865536147}
                        max_values = {'nFRing': 3, 'VR1_A': 2438706322851535.0, 'n6ARing': 6, 'SMR_VSA7': 167.61430607777726, 'nBase': 5, 'BCUTs-1h': 8.010515446905131,
                        'BCUTdv-1l': 2.8123468888969216, 'BCUTZ-1l': 5.809984371085188, 'GATS5i': 1.4518815108596792, 'BCUTdv-1h': 7.017362493629758}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)
        agree_draw_smiles = st.checkbox("Draw chemical structure")
        if agree_draw_smiles:
            smile_code = st_ketcher()
            st.markdown(f"SMILES: {smile_code}")
        st.write("---")
        st.subheader("5-HT2A in pharmacology")
        st.write("The 5-HT2A receptors are abundant in cortical and limbic areas, the basal ganglia, and mesolimbic pathways. They play a crucial role in various physiological and behavioral responses.")
        st.write("Strong evidence from selective pharmacological and genetic studies supports the 5-HT2A receptor's involvement in animal responses, including head twitch responses, drug discrimination signals, hyperthermia, and changes in exploratory behavior. Moreover, there is substantial evidence linking the 5-HT2A receptor to the hallucinogenic effects of psychedelic drugs in humans. The head twitch response induced by a 5-HT2A receptor agonist has predictive value and involves the phosphoinositide signaling pathway, primarily in the prefrontal cortex of the forebrain.")
        st.write("Many newer antipsychotic drugs, such as aripiprazole and brexpiprazole, exhibit high affinity for the 5-HT2A receptor, though their action may also involve other receptors, like the 5-HT1A receptor.")
        col1, col2 = st.columns(2)
        with col1:
            image = Image.open('images_app/aripiprazol_app.png')
            st.image(image, caption='Aripiprazole', use_column_width=True)
        with col2:
            image = Image.open('images_app/breksipiprazol_app.png')
            st.image(image, caption='Brexpiprazole', use_column_width=True)
        st.write("Recent preclinical findings indicate that 5-HT2A receptor antagonists enhance the antidepressant effects of 5-HT2A uptake inhibitors and increase extracellular serotonin levels. This suggests potential adjunctive use of 5-HT2A receptor antagonists in treatment-resistant depression. The exact mechanism underlying this effect is not fully understood but may involve inhibitory feedback mediated by the 5-HT2A receptor on serotonin neurons.")
        st.write("Interestingly, the 5-HT2A receptor is associated with impulse control. Selective 5-HT2A receptor antagonists consistently reduce impulsivity in animal models and have been confirmed by psycho-pharmacological studies in humans, such as the use of quetiapine, a non-selective 5-HT2A receptor antagonist.")
        st.write("These findings raise the possibility that 5-HT2A receptor antagonists/inverse agonists, along with novel agents like ebselen, could be useful in controlling impulse control disorders.")
        

    elif sub_option == "5-HT2B receptor":
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            predictions = best_model_5HT2B.predict(descriptors_value_df)
                        st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                        prediction_float = round(float(predictions), 3)
                        st.write("pKi value for 5-HT2B serotonin receptor: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['VSA_EState4', 'AATS6i', 'IC2', 'ATSC6s', 'AATSC6p', 
                                                        'AMID_O', 'GATS6c', 'AATS7se', 'MDEC-23', 'MATS7se']
                        min_values = {'VSA_EState4': -2.751429417977417, 'AATS6i': 149.59935464493668, 'IC2': 2.4806821149663847, 'ATSC6s': -112.09948979591836,
                        'AATSC6p': -0.2129412957609489, 'AMID_O': 0.0, 'GATS6c': 0.1545832295965297, 'AATS7se': 6.718464, 'MDEC-23': 0.0, 'MATS7se': -1.8934399907692108}
                        max_values = {'VSA_EState4': 21.24730570625863, 'AATS6i': 187.08515599534223, 'IC2': 5.356371641640362, 'ATSC6s': 468.17898022892814,
                        'AATSC6p': 0.2242242923218947, 'AMID_O': 0.7540191188648945, 'GATS6c': 1.7990059889627592, 'AATS7se': 9.9144,
                        'MDEC-23': 51.501855527054175, 'MATS7se': 0.9962877708402362}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)
        agree_draw_smiles = st.checkbox("Draw chemical structure")
        if agree_draw_smiles:
            smile_code = st_ketcher()
            st.markdown(f"SMILES: {smile_code}")
        st.write("---")
        st.subheader("5-HT2B in pharmacology")
        st.write("The 5-HT2B receptor, a subtype of serotonin receptors, plays a significant role in regulating various physiological and behavioral processes. Interestingly, studies using knockout mice lacking the 5-HT2B receptor have revealed striking phenotypes, including deficits in sensory-motor gating, social interactions, attention, learning, and memory. Additionally, these mice display increased impulsivity and altered sleep patterns. Some of these phenotypic changes can be mimicked by blocking the 5-HT2B receptor using a selective antagonist called RS127445. These findings highlight the importance of the 5-HT2B receptor in modulating complex behaviors and cognitive functions, such as sensory processing, social behavior, attention, learning, memory, impulsivity, and sleep. Further research is needed to fully understand the underlying mechanisms and potential therapeutic implications of targeting the 5-HT2B receptor in these domains.")
        col1, col2, col3 = st.columns([1,2,1])
        with col2:
            image = Image.open('images_app/5-HT2B.png')
            st.image(image, use_column_width=True)

    
    elif sub_option == "5-HT2C receptor":
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            predictions = best_model_5HT2C.predict(descriptors_value_df)
                        st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                        prediction_float = round(float(predictions), 3)
                        st.write("pKi value for 5-HT2C serotonin receptor: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['SlogP_VSA2', 'SlogP_VSA1', 'Diameter', 'SMR_VSA6', 'nFAHRing', 
                                                        'nAHRing', 'nFaHRing', 'MDEO-11', 'MDEC-33', 'JGI4']
                        min_values = {'SlogP_VSA2': 0.0, 'SlogP_VSA1': 0.0, 'Diameter': 5, 'SMR_VSA6': 0.0, 'nFAHRing': 0,
                        'nAHRing': 0, 'nFaHRing': 0, 'MDEO-11': 3.4657242157757293e-06, 'MDEC-33': 0.0, 'JGI4': 0.0062222222222222}
                        max_values = {'SlogP_VSA2': 517.9559660933347, 'SlogP_VSA1': 33.15804067660034, 'Diameter': 100000000, 'SMR_VSA6': 112.31699339671962,
                        'nFAHRing': 2, 'nAHRing': 19, 'nFaHRing': 2, 'MDEO-11': 18.414657040197348, 'MDEC-33': 33.010349012228495, 'JGI4': 0.0826666666666666}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)
        agree_draw_smiles = st.checkbox("Draw chemical structure")
        if agree_draw_smiles:
            smile_code = st_ketcher()
            st.markdown(f"SMILES: {smile_code}")
        st.write("---")
        st.subheader("5-HT2C in pharmacology")
        st.write("The expression of the 5-HT2C receptor is mainly confined to the central nervous system. Abundant distribution of 5-HT2C receptors occurs in cortical and limbic areas, the basal ganglia, mesolimbic pathways, and the hypothalamus.")
        st.write("Early neuropharmacological studies in rodents using non-selective agents associated activation of the 5-HT2C receptor with various behavioral and physiological effects, such as reduced locomotion, decreased food intake, anxiety, and increased body temperature. With the introduction of more selective agonists (e.g., CP-809101, lorcaserin) and antagonists (e.g., SB242084), along with studies using 5-HT2C knockout mice, the range of functions attributed to the 5-HT2C receptor expanded to include compulsive drug-seeking and feeding behavior, central control of energy homeostasis, orofacial dyskinesia, vigilance, and seizure threshold control.")
        col1, col2 = st.columns(2)
        with col1:
            st.write("A selective 5-HT2C receptor agonist, lorcaserin, has undergone clinical trials and been approved for the treatment of obesity (in combination with lifestyle modifications). Although the overall magnitude of lorcaserin's effect appears modest, individual differences in response indicate that patient stratification could be utilized to optimize treatment outcomes in the future.")
            st.write("The 5-HT2C receptor, with its distinct expression pattern and involvement in various physiological and behavioral processes, presents a promising target for therapeutic interventions in conditions such as obesity. Further research and understanding of the receptor's complex functions will pave the way for potential advancements in the treatment of related disorders.")
        with col2:
            image = Image.open('images_app/5-HT2C.png')
            st.image(image, use_column_width=True)            
        
        
    elif sub_option == "5-HT3 receptor":
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            predictions = best_model_5HT3.predict(descriptors_value_df)
                        st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                        prediction_float = round(float(predictions), 3)
                        st.write("pKi value for 5-HT3 serotonin receptor: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['AATSC8d', 'nFRing', 'nBridgehead', 'NaasN', 'AATS7p', 'AATSC6Z',
                                                        'GATS8s', 'AATS6d', 'MINaaN', 'MAXaaN']
                        min_values = {'AATSC8d': -0.2262323943661972, 'nFRing': 0, 'nBridgehead': 0, 'NaasN': 0, 'AATS7p': 0.6767748409994595,
                        'AATSC6Z': -16.57908163265306, 'GATS8s': 0.0, 'AATS6d': 1.40625, 'MINaaN': 3.5102435279667423, 'MAXaaN': 3.5102435279667423}
                        max_values = {'AATSC8d': 0.5625, 'nFRing': 3, 'nBridgehead': 4, 'NaasN': 2, 'AATS7p': 1.8218280345317268,
                        'AATSC6Z': 9.779639889196677, 'GATS8s': 3.378510673341657, 'AATS6d': 3.761904761904762, 'MINaaN': 5.141287792894936, 'MAXaaN': 5.141287792894936}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)
        agree_draw_smiles = st.checkbox("Draw chemical structure")
        if agree_draw_smiles:
            smile_code = st_ketcher()
            st.markdown(f"SMILES: {smile_code}")                
        st.write("---")
        st.subheader("5-HT3 in pharmacology")
        col1, col2 = st.columns([1,3])
        with col1:
            image = Image.open('images_app/5-HT3.png')
            st.image(image, use_column_width=True)              
        with col2:
            st.write("The 5-HT3 receptors are abundantly found in cortical and limbic areas, as well as in the dorsal vagal complex of the vagus nerve.")
            st.write("In addition to their well-established role in controlling nausea and vomiting (involving both central and peripheral mechanisms), the 5-HT3 receptors have been associated with various behavioral effects. These range from changes in anxiety and cognitive function to altered pain processing and sensitivity to addictive substances.")
        st.write("Recent evidence has linked 5-HT3 receptors to emotions, behaviors, and cognitive functions, particularly through findings related to vortioxetine. Vortioxetine exhibits atypical pharmacology at the 5-HT3 receptor, initially acting as an agonist followed by long-lasting, irreversible receptor inhibition. These effects of vortioxetine are associated with its high affinity for binding sites of 5-HT1A, 5-HT1B, 5-HT1D, 5-HT7, and serotonin transporters. Vortioxetine has significant antidepressant and pro-cognitive activity in both rodent models and clinical studies and is currently marketed as an antidepressant with cognitive-enhancing properties. Despite its polymodal pharmacology, blocking 5-HT3 receptors is believed to play a significant role in the mechanism of action of vortioxetine, particularly regarding cognitive function.")
        st.write("The 5-HT3 receptors' involvement in various physiological and behavioral processes, including nausea and vomiting control, anxiety, cognitive function, pain processing, and sensitivity to addictive substances, highlights their significance as potential targets for therapeutic interventions. Further research and understanding of the complex functions of 5-HT3 receptors hold promise for the development of novel treatments addressing these conditions.")
        
    elif sub_option == "5-HT4 receptor":
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            predictions = best_model_5HT4.predict(descriptors_value_df)
                        st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                        prediction_float = round(float(predictions), 3)
                        st.write("pKi value for 5-HT4 serotonin receptor: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['C1SP3', 'SMR_VSA4', 'NssCH2', 'ATSC0c', 'C3SP3', 'ATSC8c', 'AATS3i', 'ATSC1dv', 'AATS4i', 'ATSC4c']
                        min_values ={'C1SP3': 0, 'SMR_VSA4': 0.0,  'NssCH2': 0,  'ATSC0c': 0.2838127021284238,
                         'C3SP3': 0, 'ATSC8c': -1.152247769776361, 'AATS3i': 148.91780824944055, 'ATSC1dv': 4.117167133670763,
                         'AATS4i': 150.77796634225993, 'ATSC4c': -0.6393548403766612}
                        max_values = {'C1SP3': 20, 'SMR_VSA4': 41.42534232312975, 'NssCH2': 38, 'ATSC0c': 3.818084087512284,
                         'C3SP3': 7, 'ATSC8c': 0.9405402774301158, 'AATS3i': 166.23324570851813, 'ATSC1dv': 238.99467400128788,
                         'AATS4i': 177.48398206484842, 'ATSC4c': 1.6028283796421825}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)
        agree_draw_smiles = st.checkbox("Draw chemical structure")
        if agree_draw_smiles:
            smile_code = st_ketcher()
            st.markdown(f"SMILES: {smile_code}")                
        st.write("---")
        st.subheader("5-HT4 in pharmacology")
        st.write("The 5-HT4 receptor primarily resides in the brain, particularly in the hippocampus. It plays a crucial role in cognition, with 5-HT4 receptor agonists demonstrating pro-cognitive effects in various species and memory paradigms.")            
        st.write("In addition to their well-established role in controlling nausea and vomiting (involving both central and peripheral mechanisms), the 5-HT3 receptors have been associated with various behavioral effects. These range from changes in anxiety and cognitive function to altered pain processing and sensitivity to addictive substances.")
        st.write("Additionally, these receptors are implicated in regulating feeding behavior, as observed in studies with knockout mice.")
        st.write("Despite promising cognitive benefits, the clinical development of 5-HT4 receptor agonists may face challenges, as evidenced by potential adverse effects and the need for further studies to assess general tolerability.")
        



    elif sub_option == "5-HT5A receptor":
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            predictions = best_model_5HT5A.predict(descriptors_value_df)
                        st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                        prediction_float = round(float(predictions), 3)
                        st.write("pKi value for 5-HT5A serotonin receptor: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['nBase', 'n10FRing', 'SpDiam_A', 'BCUTv-1h', 'BCUTZ-1l', 'ATSC3s', 'BCUTc-1l', 
                                                        'BCUTi-1l', 'MINdssC', 'GGI4']
                        min_values = {'nBase': 0, 'n10FRing': 0, 'SpDiam_A': 4.417866110045354, 'BCUTv-1h': 20.860460092488022,
                        'BCUTZ-1l': 0.9918795477066924, 'ATSC3s': -64.16436554898094, 'BCUTc-1l': -0.5362676824244084,
                        'BCUTi-1l': 10.3110261839478, 'MINdssC': -2.756944444444444, 'GGI4': 0.3772222222222222}
                        max_values = {'nBase': 5, 'n10FRing': 2, 'SpDiam_A': 6.105734449269413, 'BCUTv-1h': 32.51604863308039,
                        'BCUTZ-1l': 5.789018651182438, 'ATSC3s': 377.3572708000317, 'BCUTc-1l': -0.3025224914355429,
                        'BCUTi-1l': 11.06952158145765, 'MINdssC': 1.17587962962963, 'GGI4': 4.712222222222223}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)
        agree_draw_smiles = st.checkbox("Draw chemical structure")
        if agree_draw_smiles:
            smile_code = st_ketcher()
            st.markdown(f"SMILES: {smile_code}")                
        st.write("---")
        st.subheader("5-HT5A in pharmacology")
        col1, col2 = st.columns([3,1])
        with col1:
            st.write("Both subtypes of the 5-HT5 receptor, 5-HT5A and 5-HT5B, have been detected in the central nervous system of various species. However, the full-length 5-HT5B receptor is not expressed in humans, which diminishes its therapeutic relevance as a target. Currently, limited information is available regarding the functions of the central nervous system 5-HT5A receptor due to the lack of appropriately selective pharmacological tools. Several drugs with 5-HT5A receptor antagonist activity have been studied in rodents for their potential anxiolytic, antidepressant, and antipsychotic effects. However, a consistent behavioral effect of blocking the 5-HT5A receptor has not been identified thus far, although a possible relationship with memory has been observed.")
        st.write("The 5-HT5A receptor remains an area of ongoing research, and further studies are needed to uncover its precise functions and potential therapeutic implications.")
        with col2:
            image = Image.open('images_app/5-HT5A.png')
            st.image(image, use_column_width=True)   
        

    elif sub_option == "5-HT6 receptor":
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            predictions = best_model_5HT6.predict(descriptors_value_df)
                        st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                        prediction_float = round(float(predictions), 3)
                        st.write("pKi value for 5-HT6 serotonin receptor: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['NddssS', 'BCUTd-1h', 'PEOE_VSA2','PEOE_VSA12', 'PEOE_VSA1', 
                                 'n9FaRing', 'VR1_A', 'MINssCH2', 'JGI5', 'MINsssN']
                        min_values = {'NddssS': 0, 'BCUTd-1h': 3.077886962639631, 'PEOE_VSA2': 0.0, 'PEOE_VSA12': 0.0, 'PEOE_VSA1': 0.0, 'n9FaRing': 0,
                        'VR1_A': 36.7685702901585, 'MINssCH2': -1.4125000468962576, 'JGI5': 0.0, 'MINsssN': 0.5848844954648524}
                        max_values = {'NddssS': 3, 'BCUTd-1h': 4.123286793361294, 'PEOE_VSA2': 32.59097764862998, 'PEOE_VSA12': 26.638119238973665,
                        'PEOE_VSA1': 31.40487476392844, 'n9FaRing': 2, 'VR1_A': 2353315.256893448, 'MINssCH2': 1.2974889019628308, 'JGI5': 0.0457765151515151,
                        'MINsssN': 2.81188171747654}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)
        agree_draw_smiles = st.checkbox("Draw chemical structure")
        if agree_draw_smiles:
            smile_code = st_ketcher()
            st.markdown(f"SMILES: {smile_code}")                
        st.write("---")
        st.subheader("5-HT6 in pharmacology")
        col1, col2 = st.columns([3,1])
        with col1:
            st.write("Expression of the 5-HT6 receptor is limited to the central nervous system. These receptors are abundantly present in cortical and limbic regions, as well as in the basal ganglia.")
            st.write("Although the initial ligands for the 5-HT6 receptor lacked selectivity and brain penetration, these challenges have been largely overcome by recently developed compounds, including the agonists WAY181187 and WAY208466, as well as the antagonists LuAE58054 (idalopirdine), SB258585, and SB399885. The use of such compounds in behavioral models of rats has allowed for the identification of the role of the 5-HT6 receptor in a range of central nervous system functions, including learning and memory, feeding, addictive behaviors, and seizure control.")
            st.write("In particular, there is consistent evidence that 5-HT6 receptor antagonists improve learning and memory in preclinical models, ranging from novel object recognition and social recognition to spatial memory tasks.")
            st.write("Partial or inverse agonists of the 5-HT6 receptor represent potential candidates for future research on Alzheimer's disease and other disorders characterized by cognitive deficits.")
        with col2:
            image = Image.open('images_app/5-HT6.png')
            st.image(image, use_column_width=True)  

    
    elif sub_option == "5-HT7 receptor":
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            predictions = best_model_5HT7.predict(descriptors_value_df)
                        st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                        prediction_float = round(float(predictions), 3)
                        st.write("pKi value for 5-HT7 serotonin receptor: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['BCUTd-1h', 'GATS8d', 'PEOE_VSA9', 'AATS6v', 'AATS7v', 'MAXdO', 
                                                        'AATS6d', 'GATS8p', 'AATSC4d', 'Xch-6d']
                        min_values = {'BCUTd-1h': 3.0379937736676954, 'GATS8d': 0.0, 'PEOE_VSA9': 0.0, 'AATS6v': 72.93864269488313,
                        'AATS7v': 31.083744430930263, 'MAXdO': 9.099537037037038, 'AATS6d': 1.4137931034482758, 'GATS8p': 0.0,
                        'AATSC4d': -0.3652949245541836, 'Xch-6d': 0.0}
                        max_values = {'BCUTd-1h': 4.128970079454183, 'GATS8d': 1.862385321100917, 'PEOE_VSA9': 43.04204896991582, 'AATS6v': 280.56396703468283,
                        'AATS7v': 256.1013420435136, 'MAXdO': 14.188116512215537, 'AATS6d': 3.761904761904762, 'GATS8p': 4.324183757840702, 'AATSC4d': 0.1093325407608695,
                        'Xch-6d': 4.768941808213231}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")                      
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)  
        agree_draw_smiles = st.checkbox("Draw chemical structure")
        if agree_draw_smiles:
            smile_code = st_ketcher()
            st.markdown(f"SMILES: {smile_code}")                
        st.write("---")
        st.subheader("5-HT7 in pharmacology")
        st.write("The 5-HT7 receptors are abundantly present in the suprachiasmatic nucleus.")
        st.write("Currently, selective agonists (AS-19, LP-211, E-55888) as well as antagonists (SB-258719, SB269970) are available. Previous studies that identified 8-OH-DPAT-induced hypothermia suggested the involvement of the 5-HT1A receptor before this drug was also recognized as having agonistic properties at the 5-HT7 receptor. It is now evident that 5-HT7 receptor agonists, such as AS-19, induce hypothermia in various species, and both 5-HT1A and 5-HT7 receptors are involved in the action of 8-OH-DPAT.")
        col1, col2 = st.columns([1,2])
        with col1:
            image = Image.open('images_app/5-HT7.png')
            st.image(image, use_column_width=True)
        with col2:
            st.write("Evidence suggests that the 5-HT7 receptor is not only associated with temperature control but also plays a role in a wide range of central nervous system processes, including the regulation of circadian rhythms, mood, learning and memory, seizure threshold, pain processing, and addiction mechanisms.")
            st.write("Other antidepressant and antipsychotic drugs (tricyclics, lurasidone, aripiprazole, etc.) used in clinical practice possess antagonist properties at the 5-HT7 receptor.")
                
     
elif selected == "SERT":
    calc = Calculator(descriptors, ignore_3D=True)
    calc.descriptors = [d for d in calc.descriptors if str(d) in descriptors_for_QSAR] 
    st.markdown('<h1 class="text-second-title">Predictions for the serotonin transporter</h1>', unsafe_allow_html=True)
    smiles_input = st.text_input("Input SMILES", key="text")
    col1, col2 = st.columns(2)
    if smiles_input:
        try:
            molecule = Chem.MolFromSmiles(smiles_input)
            if molecule:
                img = Draw.MolToImage(molecule)
                with col1:
                    st.image(img, caption='Chemical structure', use_column_width=True)
            else:
                pass
        except Exception as e:
            st.error(f"Wystąpił błąd: {str(e)}")
    if smiles_input:
        try:
            molecule = Chem.MolFromSmiles(smiles_input)
            if molecule is not None:
                descriptors_value = calc.pandas([molecule])
                descriptors_value_df = pd.DataFrame(descriptors_value)
                for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                    descriptors_value_df[column] = 0
                with col2:
                    with st.spinner('Calculation in progress'):
                        predictions = best_model_sert.predict(descriptors_value_df)
                    st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                    prediction_float = round(float(predictions), 3)
                    st.write("pKi value for serotonin transporter: ", f'<span style="color: #5d93a3;">{ prediction_float}</span>', unsafe_allow_html=True)
                    st.button("Clear SMILES", on_click=clear_text)
                    list_of_important_descriptors = ['SIC2', 'GATS5c', 'IC2', 'nBase', 'ATSC2d', 'CIC1', 
                                 'SLogP', 'n10FaRing', 'SlogP_VSA1', 'PEOE_VSA3']
                    min_values = {'SIC2': 0.4663428534035417, 'GATS5c': 0.1522646188176759, 'IC2': 2.754636215098623, 'nBase': 0,
                    'ATSC2d': -17.029333333333327, 'CIC1': 0.7393028412041107, 'SLogP': -1.0360999999999996,
                    'n10FaRing': 0, 'SlogP_VSA1': 0.0, 'PEOE_VSA3': 0.0}
                    max_values = {'SIC2': 0.9484501704693528, 'GATS5c': 1.9718023087752352, 'IC2': 5.48361752790408, 'nBase': 3, 'ATSC2d': 19.849172805216764,
                    'CIC1': 4.040885762787749, 'SLogP': 11.12039999999999, 'n10FaRing': 2, 'SlogP_VSA1': 48.53088597686407, 'PEOE_VSA3': 30.5236297402646}
                    normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                    values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                    values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                    values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                    labels = normalized_descriptors_df[list_of_important_descriptors].columns
                    desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                    values_1 += values_1[:1]
                    values_2 += values_2[:1]
                    values_2 = [-1 if value < -1 else value for value in values_2]
                    values_2 = [1.5 if value > 1.5 else value for value in values_2]
                    values_3 += values_3[:1]
                    num_labels = len(labels)
                    angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                    angles += angles[:1]
                    fig = plt.figure(figsize=(8,8))
                    ax = fig.add_subplot(111, polar=True)
                    color_1 = '#A6A6A6'
                    color_2 = '#4282AA'
                    ax.plot(angles, values_1, color=color_1, label="training set")
                    ax.fill(angles, values_1, alpha=0.25, color=color_1)
                    ax.plot(angles, values_3, color="white")
                    ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                    ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                    ax.fill(angles, values_2, alpha=0)
                    ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                    ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                    ax.legend()
                    plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                    plt.tight_layout()
                    st.pyplot(fig)    
                    if desc_condition > 6:
                        st.write("Compound is under applicability domain")
                    else:
                        st.write("Compound is not under applicability domain, prediction may be inaccurate")
            else:
                st.write("Invalid SMILES")
        except Exception as e:
            st.write("Error:", e)
    agree_draw_smiles = st.checkbox("Draw chemical structure")
    if agree_draw_smiles:
        smile_code = st_ketcher()
        st.markdown(f"SMILES: {smile_code}")
    st.write("---")
    st.subheader("SERT in pharmacology")
    st.write("The serotonin transporter is a protein responsible for transporting serotonin, a neurotransmitter, between the synapses of neurons. It is a key component in the regulation of serotonin levels in the brain and is involved in numerous functions, such as the regulation of mood, appetite, sleep and cognitive processes. The serotonin transporter works by taking serotonin from the synaptic space and transporting it back to the presynaptic neuron, where it can be re-released or degraded by enzymes. Malfunctions in the serotonin transporter can contribute to various mental illnesses, such as depression and anxiety disorders.")
    col1, col2 = st.columns([3,1])
    with col1:
        st.write("The serotonin transporter is also a major target of SSRI (selective serotonin reuptake inhibitor) drugs, which increase the concentration of serotonin in the synaptic space, improving communication between neurons. By blocking the serotonin transporter, SSRI drugs can increase the availability of serotonin, which can help improve mood and reduce symptoms of depression and other psychiatric disorders.")
    with col2:
        image = Image.open('images_app/SERT.png')
        st.image(image, use_column_width=True)        
    
            
elif selected == "HIA":
    st.markdown('<h1 class="text-second-title">Class prediction of human instestinal absorption (HIA)</h1>', unsafe_allow_html=True)
    descriptors_for_QSPR = pd.read_table("descriptors_QSPR.txt", header = None)[0][0]
    calc = Calculator(descriptors, ignore_3D=True)
    calc.descriptors = [d for d in calc.descriptors if str(d) in descriptors_for_QSPR] 
    with st.container():
        st.write("Prediction based on single SMILES")
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    calc = Calculator(descriptors, ignore_3D=True)
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            probability = best_model_hia_classification.predict_proba(descriptors_value_df)
                            if probability[0][1] >= 0.775:
                                st.write("Good intestinal permeability")
                            else:
                                regression_predictions = best_model_hia_regression.predict(descriptors_value_df)
                                if regression_predictions >= 90:
                                    st.write("Good intestinal permeability")
                                else:
                                    st.write("Bad intestinal permeability")
                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['MATS1dv', 'PEOE_VSA7', 'MATS4s', 'RPCG', 'AATSC4p', 'AMID_X', 
                                        'ATSC6i', 'AATSC3d', 'MID_N', 'GATS4i']
                        min_values = {'MATS1dv': -0.1501196994477234, 'PEOE_VSA7': 6.923737199690624, 'MATS4s': -0.2164223259483269,
                        'RPCG': 0.0532270790025828, 'AATSC4p': -0.1272635231618815, 'AMID_X': 0.0, 'ATSC6i': -64.61742537264104,
                        'AATSC3d': -0.2011699790527768, 'MID_N': 0.0, 'GATS4i': 0.6972267893534747}
                        max_values = {'MATS1dv': 0.5060337892196302, 'PEOE_VSA7': 81.11421727488754, 'MATS4s': 0.3658034291146026,
                        'RPCG': 0.3065626346663486, 'AATSC4p': 0.0931950415924, 'AMID_X': 0.282628146929378,
                        'ATSC6i': 45.86523543779616, 'AATSC3d': 0.1227586206896551, 'MID_N': 14.58274965700862,
                        'GATS4i': 1.358378838745984}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)  
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)        
    agree_draw_smiles = st.checkbox("Draw chemical structure")
    if agree_draw_smiles:
        smile_code = st_ketcher()
        st.markdown(f"SMILES: {smile_code}")      
    with st.container():
        st.write("---")
        st.write('Predictions of HIA based on uploaded file')
        uploaded_file = st.file_uploader('CSV file')
        if uploaded_file is not None:
            calc = Calculator(descriptors, ignore_3D=True)
            data_file = pd.read_csv(uploaded_file)
            st.write(data_file.head())
            descriptor = []
            smiles = []
            for smi in data_file["smiles"]:
                mols = [Chem.MolFromSmiles(smi)]
                smiles.append(smi)
                descriptor.append((calc.pandas(mols)))
            descriptor_df = pd.DataFrame(np.vstack(descriptor))
            descriptor_df.columns = descriptor[0].columns
            descriptor_df = descriptor_df.astype(float)
            descriptor_df = descriptor_df.fillna(0)
            list_of_important_descriptors = ['MATS1dv', 'PEOE_VSA7', 'MATS4s', 'RPCG', 'AATSC4p', 'AMID_X', 
                                        'ATSC6i', 'AATSC3d', 'MID_N', 'GATS4i']
            min_values = {'MATS1dv': -0.1501196994477234, 'PEOE_VSA7': 6.923737199690624, 'MATS4s': -0.2164223259483269,
                        'RPCG': 0.0532270790025828, 'AATSC4p': -0.1272635231618815, 'AMID_X': 0.0, 'ATSC6i': -64.61742537264104,
                        'AATSC3d': -0.2011699790527768, 'MID_N': 0.0, 'GATS4i': 0.6972267893534747}
            max_values = {'MATS1dv': 0.5060337892196302, 'PEOE_VSA7': 81.11421727488754, 'MATS4s': 0.3658034291146026,
                        'RPCG': 0.3065626346663486, 'AATSC4p': 0.0931950415924, 'AMID_X': 0.282628146929378,
                        'ATSC6i': 45.86523543779616, 'AATSC3d': 0.1227586206896551, 'MID_N': 14.58274965700862,
                        'GATS4i': 1.358378838745984}
            normalized_descriptors_df = (descriptor_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
            data = pd.DataFrame()
            if st.button("Calculate Predictions"):
                data = pd.DataFrame(columns = ["HIA_class"])
                for i in range(len(descriptor_df)):
                    compound = pd.DataFrame(descriptor_df.iloc[i]).transpose()
                    probability = best_model_hia_classification.predict_proba(compound)
                    if probability[0][1] >= 0.775:
                        data.loc[i, "HIA_class"] = "Good intestinal permeability"
                    else:
                        regression_predictions = best_model_hia_regression.predict(compound)
                        if regression_predictions >= 90:
                            data.loc[i, "HIA_class"] = "Good intestinal permeability"
                        else:
                            data.loc[i, "HIA_class"] = "Poor intestinal permeability"
                    values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                    values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                    values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                    desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                    if desc_condition > 6:
                        data.loc[i, "applicability_domain"] = True  
                    else:
                        data.loc[i, "applicability_domain"] = False
                data_new = pd.concat([data_file, data], axis = 1)
                if data_file.empty:
                    st.write("No predictions to download.")
                else:
                    st.download_button(label="Download predictions as csv file", data=data_new.to_csv(index=False), file_name="predictions.csv")                        
    st.write("---")
    st.subheader("Human intestinal absorption")
    st.write( "HIA, or Human Intestinal Absorption, is the prediction of the ability of a drug or compound to be absorbed into the bloodstream after oral administration. It is a crucial factor in drug development as the oral route is the most preferred and natural way of drug dosing. By accurately predicting HIA, researchers can screen and prioritize potential oral drug candidates early in the drug discovery process, saving time and resources. AI-based approaches, such as the one described in this study focusing on serotonergic molecules, offer a promising method for predicting HIA and improving the efficiency of drug development.")
    image = Image.open('images_app/HIA.png')
    st.image(image, use_column_width=True) 
    st.write("The AI-based system was developed with the involvement of Klaudia Klimończyk as part of her master's degree project.")
    
    
elif selected == "BBB":
    st.markdown('<h1 class="text-second-title">Class prediction of Blood-Brain Barrier penetration (BBB)</h1>', unsafe_allow_html=True)
    descriptors_for_QSPR_BBB = pd.read_table("descriptors_QSPR_BBB.txt", header = None)[0][0]
    calc = Calculator(descriptors, ignore_3D=True)
    calc.descriptors = [d for d in calc.descriptors if str(d) in descriptors_for_QSPR_BBB] 
    with st.container():
        st.write("Prediction based on single SMILES")
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    calc = Calculator(descriptors, ignore_3D=True)
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            probability = best_model_BBB.predict_proba(descriptors_value_df)
                            if probability[0][1] > 0.8638:
                                st.write("Good blood-brain penetration")
                            else:
                                st.write('Bad blood-brain penetation')

                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['TopoPSA(NO)', 'MDEC-33', 'SlogP_VSA10', 'MAXsOH', 'AATS1i', 'Xch-5dv', 
                   'MATS2s', 'MINdssC', 'EState_VSA4', 'AATS7d']
                        min_values = {'TopoPSA(NO)': 0.0, 'MDEC-33': 0.0, 'SlogP_VSA10': 0.0, 'MAXsOH': 0.0, 'AATS1i': 0.0, 'Xch-5dv': 0.0,
                        'MATS2s': -0.7, 'MINdssC': -2.942453703703704, 'EState_VSA4': 0.0, 'AATS7d': 0.0}
                        max_values = {'TopoPSA(NO)': 662.4100000000001, 'MDEC-33': 50.8420005596026, 'SlogP_VSA10': 52.68498057209785, 'MAXsOH': 13.713664854183364,
                        'AATS1i': 211.24006281, 'Xch-5dv': 1.9072548561523543, 'MATS2s': 0.7938114236332964, 'MINdssC': 1.5624537037037038,
                        'EState_VSA4': 123.95698427278327, 'AATS7d': 3.5225225225225225}                        
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig) 
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)        
    agree_draw_smiles = st.checkbox("Draw chemical structure")
    if agree_draw_smiles:
        smile_code = st_ketcher()
        st.markdown(f"SMILES: {smile_code}")      
    with st.container():
        st.write("---")
        st.write('Predictions of BBB penetration based on uploaded file')
        uploaded_file = st.file_uploader('CSV file')
        if uploaded_file is not None:
            calc = Calculator(descriptors, ignore_3D=True)
            data_file = pd.read_csv(uploaded_file)
            st.write(data_file.head())
            descriptor = []
            smiles = []
            for smi in data_file["smiles"]:
                mols = [Chem.MolFromSmiles(smi)]
                smiles.append(smi)
                descriptor.append((calc.pandas(mols)))
            descriptor_df = pd.DataFrame(np.vstack(descriptor))
            descriptor_df.columns = descriptor[0].columns
            descriptor_df = descriptor_df.astype(float)
            descriptor_df = descriptor_df.fillna(0)
            list_of_important_descriptors = ['TopoPSA(NO)', 'MDEC-33', 'SlogP_VSA10', 'MAXsOH', 'AATS1i', 'Xch-5dv', 
                   'MATS2s', 'MINdssC', 'EState_VSA4', 'AATS7d']
            min_values = {'TopoPSA(NO)': 0.0, 'MDEC-33': 0.0, 'SlogP_VSA10': 0.0, 'MAXsOH': 0.0, 'AATS1i': 0.0, 'Xch-5dv': 0.0,
                        'MATS2s': -0.7, 'MINdssC': -2.942453703703704, 'EState_VSA4': 0.0, 'AATS7d': 0.0}

            max_values = {'TopoPSA(NO)': 662.4100000000001, 'MDEC-33': 50.8420005596026, 'SlogP_VSA10': 52.68498057209785, 'MAXsOH': 13.713664854183364,
                        'AATS1i': 211.24006281, 'Xch-5dv': 1.9072548561523543, 'MATS2s': 0.7938114236332964, 'MINdssC': 1.5624537037037038,
                        'EState_VSA4': 123.95698427278327, 'AATS7d': 3.5225225225225225} 
            normalized_descriptors_df = (descriptor_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
            data = pd.DataFrame()
            if st.button("Calculate Predictions"):
                data = pd.DataFrame(columns = ["BBB_class"])
                for i in range(len(descriptor_df)):
                    compound = pd.DataFrame(descriptor_df.iloc[i]).transpose()
                    probability = best_model_BBB.predict_proba(compound)
                    if probability[0][1] > 0.8638:
                        data.loc[i, "BBB_class"] = "Good BBB penetration"
                    else:
                        data.loc[i, "BBB_class"] = "Poor BBB penetration"
                    values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                    values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                    values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                    desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                    if desc_condition > 6:
                        data.loc[i, "applicability_domain"] = True  
                    else:
                        data.loc[i, "applicability_domain"] = False
                data_new = pd.concat([data_file, data], axis = 1)
                if data_file.empty:
                    st.write("No predictions to download.")
                else:
                    st.download_button(label="Download predictions as csv file", data=data_new.to_csv(index=False), file_name="predictions.csv")                        
    st.write("---")
    st.subheader("Blood-Brain Barrier")
    st.write( "The blood-brain barrier (BBB) is a specialized structure of blood vessels present within the central nervous system. It acts as a regulatory filter, controlling the flow of substances between the blood and nervous tissue. The BBB is a significant obstacle to the free movement of chemicals from the blood to the brain, thus limiting the access of some compounds to the brain environment.")
    col1, col2 = st.columns([2, 5])
    with col2:
        st.write('In the context of pharmacotherapy, the existence of the BBB can be both a challenge and a benefit. On the one hand, for some diseases of the central nervous system, penetration of drugs through the BBB is necessary for effective treatment.')
        st.write('Side effects related to the penetration of the blood-brain barrier are another aspect. Some drugs that are desirable in peripheral areas of the body may not be suitable for the brain due to their potential adverse effects on the nervous system. Therefore, the development of drug therapies requires a balanced approach to crossing the BBB.')
    with col1:
        image = Image.open('images_app/BBB.png')
        st.image(image, use_column_width=True) 

elif selected == "Batch calculation":
    calc = Calculator(descriptors, ignore_3D=True)
    calc.descriptors = [d for d in calc.descriptors if str(d) in descriptors_for_QSAR] 
    st.markdown('<h1 class="text-second-title">Batch mode calculation</h1>', unsafe_allow_html=True)
    st.write('Predictions for all serotonergic targets based on uploaded file')
    uploaded_file = st.file_uploader('CSV file')
    if uploaded_file is not None:
        data_file = pd.read_csv(uploaded_file)
        st.write(data_file.head())
        descriptors = []
        smiles = []
        for smi in data_file["smiles"]:
            mols = [Chem.MolFromSmiles(smi)]
            desc = calc.pandas(mols)
            smiles.append(smi)
            descriptors.append(desc)
        descriptors_df = pd.DataFrame(np.vstack(descriptors))
        descriptors_df.columns = descriptors[0].columns
        descriptors_df = descriptors_df.astype(float)
        descriptors_df = descriptors_df.fillna(0)
        st.write("Select the biological targets whose affinity you want to predict")
        options = ["5-HT1A", "5-HT1B", "5-HT1D", "5-HT2A", "5-HT2B", "5-HT2C", "5-HT3", "5-HT4", "5-HT5A", "5-HT6", "5-HT7", "SERT"]
        all_selected = st.checkbox("All serotonergic targets")
        default_options = options if all_selected else []
        selected_options = st.multiselect("Choose options", options, default=default_options)
        if st.button("Make predictions"):
            with st.spinner('Calculation in progress'):
                data = pd.DataFrame()
                for option in selected_options:
                    if option == "5-HT1A":
                        predictions = best_model_5HT1A.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ['BCUTc-1l', 'MAXdO', 'MAXaaaC', 'PEOE_VSA9', 'SMR_VSA3', 'SssCH2', 
                                 'AATS4i', 'SpAbs_DzZ', 'AATS6dv', 'VSA_EState5']
                        min_values = {'BCUTc-1l': -0.7373620145398293, 'MAXdO': 9.099537037037038, 'MAXaaaC': -0.1340347251287001,
                        'PEOE_VSA9': 0.0, 'SMR_VSA3': 0.0, 'SssCH2': -0.4661989795918362, 'AATS4i': 147.56632501478904,
                        'SpAbs_DzZ': 42.05895519992669, 'AATS6dv': 0.1538461538461538, 'VSA_EState5': -7.181078813682171}
                        max_values = {'BCUTc-1l': -0.292146392415571, 'MAXdO': 14.629783178882205, 'MAXaaaC': 1.5381250000000002,
                        'PEOE_VSA9': 78.6625871223764, 'SMR_VSA3': 39.8909626482546, 'SssCH2': 23.225582666887828,
                        'AATS4i': 175.1107976137481, 'SpAbs_DzZ': 1265.278990098867, 'AATS6dv': 7.298507462686567,
                        'VSA_EState5': 8.521368790538302}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False
                        
                    elif option == "5-HT1B":
                        predictions = best_model_5HT1B.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ['MINaaaC', 'SMR_VSA6', 'SlogP_VSA1', 'NsssN', 'VSA_EState4', 
                                 'SaaaC', 'ABC', 'BCUTZ-1l', 'SlogP_VSA2', 'BCUTse-1l']
                        min_values = {'MINaaaC': 0.0430998207874111, 'SMR_VSA6': 0.0, 'SlogP_VSA1': 0.0, 'NsssN': 0, 'VSA_EState4': -3.877327289138031, 'SaaaC': 0.0,
                        'ABC': 8.485281374238573, 'BCUTZ-1l': 0.9918795477066924, 'SlogP_VSA2': 4.9839785209472085, 'BCUTse-1l': 2.3910240890446155}
                        max_values = {'MINaaaC': 1.4785704837490554, 'SMR_VSA6': 115.64600422703208, 'SlogP_VSA1': 31.57463806993713, 'NsssN': 6,
                        'VSA_EState4': 18.322880995240126, 'SaaaC': 4.997746037925896, 'ABC': 48.90440618467739, 'BCUTZ-1l': 5.772390550568072, 'SlogP_VSA2': 137.2601831474361,
                        'BCUTse-1l': 2.498598365922732}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False
         
                        
                    elif option == "5-HT1D":
                        predictions = best_model_5HT1D.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ['BCUTse-1l', 'BCUTZ-1l', 'BalabanJ', 'MAXaaCH', 'n10FARing', 'MAXaasC', 'SlogP_VSA2', 'BCUTi-1l', 'JGI3', 'EState_VSA8']
                        min_values = {'BCUTse-1l': 2.39096093838768, 'BCUTZ-1l': 0.9918795477066924, 'BalabanJ': 4.133332842407474e-07, 'MAXaaCH': 1.0679894687443388,
                        'n10FARing': 0, 'MAXaasC': 0.1587137776339895, 'SlogP_VSA2': 4.9839785209472085,
                        'BCUTi-1l': 10.300242888872871, 'JGI3': 0.0166666666666666, 'EState_VSA8': 0.0}
                        max_values = {'BCUTse-1l': 2.498598365922729, 'BCUTZ-1l': 5.772390550568072, 'BalabanJ': 2.622001753479514, 'MAXaaCH': 2.510647543251294,
                        'n10FARing': 2, 'MAXaasC': 1.6150434618291762, 'SlogP_VSA2': 111.60791765215704, 'BCUTi-1l': 10.984220997218516,
                        'JGI3': 0.0852272727272727, 'EState_VSA8': 108.07472392106264}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False
                        
                    elif option == "5-HT2A":
                        predictions = best_model_5HT2A.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ["nFRing", "VR1_A", "n6ARing", "SMR_VSA7", "nBase", "BCUTs-1h", "BCUTdv-1l", "BCUTZ-1l", "GATS5i", "BCUTdv-1h"]
                        min_values = {'nFRing': 0, 'VR1_A': 30.671729502867567, 'n6ARing': 0, 'SMR_VSA7': 0.0, 'nBase': 0, 'BCUTs-1h': 2.164333480214669,
                        'BCUTdv-1l': 0.1514613373305368, 'BCUTZ-1l': 5.636495147426821, 'GATS5i': 0.375280485256885, 'BCUTdv-1h': 4.048936865536147}
                        max_values = {'nFRing': 3, 'VR1_A': 2438706322851535.0, 'n6ARing': 6, 'SMR_VSA7': 167.61430607777726, 'nBase': 5, 'BCUTs-1h': 8.010515446905131,
                        'BCUTdv-1l': 2.8123468888969216, 'BCUTZ-1l': 5.809984371085188, 'GATS5i': 1.4518815108596792, 'BCUTdv-1h': 7.017362493629758}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False
                        
                    elif option == "5-HT2B":
                        predictions = best_model_5HT2B.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ['VSA_EState4', 'AATS6i', 'IC2', 'ATSC6s', 'AATSC6p', 
                                                        'AMID_O', 'GATS6c', 'AATS7se', 'MDEC-23', 'MATS7se']
                        min_values = {'VSA_EState4': -2.751429417977417, 'AATS6i': 149.59935464493668, 'IC2': 2.4806821149663847, 'ATSC6s': -112.09948979591836,
                        'AATSC6p': -0.2129412957609489, 'AMID_O': 0.0, 'GATS6c': 0.1545832295965297, 'AATS7se': 6.718464, 'MDEC-23': 0.0, 'MATS7se': -1.8934399907692108}
                        max_values = {'VSA_EState4': 21.24730570625863, 'AATS6i': 187.08515599534223, 'IC2': 5.356371641640362, 'ATSC6s': 468.17898022892814,
                        'AATSC6p': 0.2242242923218947, 'AMID_O': 0.7540191188648945, 'GATS6c': 1.7990059889627592, 'AATS7se': 9.9144,
                        'MDEC-23': 51.501855527054175, 'MATS7se': 0.9962877708402362}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False                       
                        
                    elif option == "5-HT2C":
                        predictions = best_model_5HT2C.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ['SlogP_VSA2', 'SlogP_VSA1', 'Diameter', 'SMR_VSA6', 'nFAHRing', 
                                                        'nAHRing', 'nFaHRing', 'MDEO-11', 'MDEC-33', 'JGI4']
                        min_values = {'SlogP_VSA2': 0.0, 'SlogP_VSA1': 0.0, 'Diameter': 5, 'SMR_VSA6': 0.0, 'nFAHRing': 0,
                        'nAHRing': 0, 'nFaHRing': 0, 'MDEO-11': 3.4657242157757293e-06, 'MDEC-33': 0.0, 'JGI4': 0.0062222222222222}
                        max_values = {'SlogP_VSA2': 517.9559660933347, 'SlogP_VSA1': 33.15804067660034, 'Diameter': 100000000, 'SMR_VSA6': 112.31699339671962,
                        'nFAHRing': 2, 'nAHRing': 19, 'nFaHRing': 2, 'MDEO-11': 18.414657040197348, 'MDEC-33': 33.010349012228495, 'JGI4': 0.0826666666666666}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False                       
                        
                    elif option == "5-HT3":
                        predictions = best_model_5HT3.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ['AATSC8d', 'nFRing', 'nBridgehead', 'NaasN', 'AATS7p', 'AATSC6Z',
                                                        'GATS8s', 'AATS6d', 'MINaaN', 'MAXaaN']
                        min_values = {'AATSC8d': -0.2262323943661972, 'nFRing': 0, 'nBridgehead': 0, 'NaasN': 0, 'AATS7p': 0.6767748409994595,
                        'AATSC6Z': -16.57908163265306, 'GATS8s': 0.0, 'AATS6d': 1.40625, 'MINaaN': 3.5102435279667423, 'MAXaaN': 3.5102435279667423}
                        max_values = {'AATSC8d': 0.5625, 'nFRing': 3, 'nBridgehead': 4, 'NaasN': 2, 'AATS7p': 1.8218280345317268,
                        'AATSC6Z': 9.779639889196677, 'GATS8s': 3.378510673341657, 'AATS6d': 3.761904761904762, 'MINaaN': 5.141287792894936, 'MAXaaN': 5.141287792894936}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False  

                    elif option == "5-HT4":
                        predictions = best_model_5HT4.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ['C1SP3', 'SMR_VSA4', 'NssCH2', 'ATSC0c', 'C3SP3', 'ATSC8c', 'AATS3i', 'ATSC1dv', 'AATS4i', 'ATSC4c']
                        min_values ={'C1SP3': 0, 'SMR_VSA4': 0.0,  'NssCH2': 0,  'ATSC0c': 0.2838127021284238,
                         'C3SP3': 0, 'ATSC8c': -1.152247769776361, 'AATS3i': 148.91780824944055, 'ATSC1dv': 4.117167133670763,
                         'AATS4i': 150.77796634225993, 'ATSC4c': -0.6393548403766612} 
                        max_values = {'C1SP3': 20, 'SMR_VSA4': 41.42534232312975, 'NssCH2': 38, 'ATSC0c': 3.818084087512284,
                         'C3SP3': 7, 'ATSC8c': 0.9405402774301158, 'AATS3i': 166.23324570851813, 'ATSC1dv': 238.99467400128788,
                         'AATS4i': 177.48398206484842, 'ATSC4c': 1.6028283796421825}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False  
                        
                    elif option == "5-HT5A":
                        predictions = best_model_5HT5A.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ['nBase', 'n10FRing', 'SpDiam_A', 'BCUTv-1h', 'BCUTZ-1l', 'ATSC3s', 'BCUTc-1l', 
                                                        'BCUTi-1l', 'MINdssC', 'GGI4']
                        min_values = {'nBase': 0, 'n10FRing': 0, 'SpDiam_A': 4.417866110045354, 'BCUTv-1h': 20.860460092488022,
                        'BCUTZ-1l': 0.9918795477066924, 'ATSC3s': -64.16436554898094, 'BCUTc-1l': -0.5362676824244084,
                        'BCUTi-1l': 10.3110261839478, 'MINdssC': -2.756944444444444, 'GGI4': 0.3772222222222222}
                        max_values = {'nBase': 5, 'n10FRing': 2, 'SpDiam_A': 6.105734449269413, 'BCUTv-1h': 32.51604863308039,
                        'BCUTZ-1l': 5.789018651182438, 'ATSC3s': 377.3572708000317, 'BCUTc-1l': -0.3025224914355429,
                        'BCUTi-1l': 11.06952158145765, 'MINdssC': 1.17587962962963, 'GGI4': 4.712222222222223}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False                     
                        
                    elif option == "5-HT6":
                        predictions = best_model_5HT6.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ['NddssS', 'BCUTd-1h', 'PEOE_VSA2','PEOE_VSA12', 'PEOE_VSA1', 
                                 'n9FaRing', 'VR1_A', 'MINssCH2', 'JGI5', 'MINsssN']
                        min_values = {'NddssS': 0, 'BCUTd-1h': 3.077886962639631, 'PEOE_VSA2': 0.0, 'PEOE_VSA12': 0.0, 'PEOE_VSA1': 0.0, 'n9FaRing': 0,
                        'VR1_A': 36.7685702901585, 'MINssCH2': -1.4125000468962576, 'JGI5': 0.0, 'MINsssN': 0.5848844954648524}
                        max_values = {'NddssS': 3, 'BCUTd-1h': 4.123286793361294, 'PEOE_VSA2': 32.59097764862998, 'PEOE_VSA12': 26.638119238973665,
                        'PEOE_VSA1': 31.40487476392844, 'n9FaRing': 2, 'VR1_A': 2353315.256893448, 'MINssCH2': 1.2974889019628308, 'JGI5': 0.0457765151515151,
                        'MINsssN': 2.81188171747654}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False                      
                        
                    elif option == "5-HT7":
                        predictions = best_model_5HT7.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ['BCUTd-1h', 'GATS8d', 'PEOE_VSA9', 'AATS6v', 'AATS7v', 'MAXdO', 
                                                        'AATS6d', 'GATS8p', 'AATSC4d', 'Xch-6d']
                        min_values = {'BCUTd-1h': 3.0379937736676954, 'GATS8d': 0.0, 'PEOE_VSA9': 0.0, 'AATS6v': 72.93864269488313,
                        'AATS7v': 31.083744430930263, 'MAXdO': 9.099537037037038, 'AATS6d': 1.4137931034482758, 'GATS8p': 0.0,
                        'AATSC4d': -0.3652949245541836, 'Xch-6d': 0.0}
                        max_values = {'BCUTd-1h': 4.128970079454183, 'GATS8d': 1.862385321100917, 'PEOE_VSA9': 43.04204896991582, 'AATS6v': 280.56396703468283,
                        'AATS7v': 256.1013420435136, 'MAXdO': 14.188116512215537, 'AATS6d': 3.761904761904762, 'GATS8p': 4.324183757840702, 'AATSC4d': 0.1093325407608695,
                        'Xch-6d': 4.768941808213231}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False                      
                        
                    elif option == "SERT":
                        predictions = best_model_sert.predict(descriptors_df)
                        data[option] = predictions
                        list_of_important_descriptors = ['SIC2', 'GATS5c', 'IC2', 'nBase', 'ATSC2d', 'CIC1', 
                                    'SLogP', 'n10FaRing', 'SlogP_VSA1', 'PEOE_VSA3']
                        min_values = {'SIC2': 0.4663428534035417, 'GATS5c': 0.1522646188176759, 'IC2': 2.754636215098623, 'nBase': 0,
                        'ATSC2d': -17.029333333333327, 'CIC1': 0.7393028412041107, 'SLogP': -1.0360999999999996,
                        'n10FaRing': 0, 'SlogP_VSA1': 0.0, 'PEOE_VSA3': 0.0}
                        max_values = {'SIC2': 0.9484501704693528, 'GATS5c': 1.9718023087752352, 'IC2': 5.48361752790408, 'nBase': 3, 'ATSC2d': 19.849172805216764,
                        'CIC1': 4.040885762787749, 'SLogP': 11.12039999999999, 'n10FaRing': 2, 'SlogP_VSA1': 48.53088597686407, 'PEOE_VSA3': 30.5236297402646}
                        descriptors_value_df = descriptors_df[list_of_important_descriptors]
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        for i in range(len(descriptors_df)):
                            values_2 = normalized_descriptors_df[list_of_important_descriptors].iloc[i]
                            desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                            if desc_condition > 6:
                                data.loc[i, f"{option}_applicability_domain"] = True  
                            else:
                                data.loc[i, f"{option}_applicability_domain"] = False
                        
                data_new = pd.concat([data_file, data], axis = 1)
                if data.empty:
                    st.write("No predictions to download.")
                else:
                    st.download_button(label="Download predictions as csv file", data=data_new.to_csv(index=False), file_name="predictions.csv")   

elif selected == "Serotonergic activity":
    binary_model = AutoML(binary_model_path)
    st.markdown('<h1 class="text-second-title">Serotonergic activity predictions</h1>', unsafe_allow_html=True)
    descriptors_for_QSPR = pd.read_table("binary_active_inactive_descriptors.txt", header = None)[0][0]
    calc = Calculator(descriptors, ignore_3D=True)
    calc.descriptors = [d for d in calc.descriptors if str(d) in descriptors_for_QSPR] 
    with st.container():
        smiles_input = st.text_input("Input SMILES", key="text")
        col1, col2 = st.columns(2)
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule:
                    img = Draw.MolToImage(molecule)
                    with col1:
                        st.image(img, caption='Chemical structure', use_column_width=True)
                else:
                    pass
            except Exception as e:
                st.error(f"Wystąpił błąd: {str(e)}")
        if smiles_input:
            try:
                molecule = Chem.MolFromSmiles(smiles_input)
                if molecule is not None:
                    calc = Calculator(descriptors, ignore_3D=True)
                    descriptors_value = calc.pandas([molecule])
                    descriptors_value_df = pd.DataFrame(descriptors_value)
                    for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                        descriptors_value_df[column] = 0
                    with col2:
                        with st.spinner('Calculation in progress'):
                            probability = binary_model.predict_proba(descriptors_value_df)
                            if probability[0][1] > 0.568706:
                                st.write("Compound has a serotonergic activity")
                            else:
                                st.write('Compound may not have serotonergic activity')

                        st.button("Clear SMILES", on_click=clear_text)
                        list_of_important_descriptors = ['nBase', 'MATS1v', 'PEOE_VSA7','SlogP_VSA1', 'AXp-7dv', 'PEOE_VSA9', 'Xch-7dv','AATSC2dv', 'VSA_EState2','ATSC6v']
                        min_values = {'nBase': 0, 'MATS1v': -0.1497204202184233, 'PEOE_VSA7': 0.0, 'SlogP_VSA1': 0.0,
 'AXp-7dv': 0.0040737032751538, 'PEOE_VSA9': 0.0, 'Xch-7dv': 0.0, 'AATSC2dv': -0.7227527873894656,
 'VSA_EState2': -1.3062069767453963, 'ATSC6v': -3489.3195470956503}
                        max_values = {'nBase': 5, 'MATS1v': 0.2019023062607182, 'PEOE_VSA7': 154.4803930624426, 'SlogP_VSA1': 91.3750762459032,
 'AXp-7dv': 0.0565596697830878, 'PEOE_VSA9': 163.71203852694237, 'Xch-7dv': 10.33296316622819, 'AATSC2dv': 2.582222222222222,
 'VSA_EState2': 87.06048435105046, 'ATSC6v': 3956.085710356009}
                        normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                        values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                        values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                        values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                        labels = normalized_descriptors_df[list_of_important_descriptors].columns
                        desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                        values_1 += values_1[:1]
                        values_2 += values_2[:1]
                        values_2 = [-1 if value < -1 else value for value in values_2]
                        values_2 = [1.5 if value > 1.5 else value for value in values_2]
                        values_3 += values_3[:1]
                        num_labels = len(labels)
                        angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                        angles += angles[:1]
                        fig = plt.figure(figsize=(8,8))
                        ax = fig.add_subplot(111, polar=True)
                        color_1 = '#A6A6A6'
                        color_2 = '#4282AA'
                        ax.plot(angles, values_1, color=color_1, label="training set")
                        ax.fill(angles, values_1, alpha=0.25, color=color_1)
                        ax.plot(angles, values_3, color="white")
                        ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                        ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                        ax.fill(angles, values_2, alpha=0)
                        ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                        ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                        ax.legend()
                        plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                        plt.tight_layout()
                        st.pyplot(fig)    
                        if desc_condition > 6:
                            st.write("Compound is under applicability domain")
                        else:
                            st.write("Compound is not under applicability domain, prediction may be inaccurate")
                else:
                    st.write("Invalid SMILES")
            except Exception as e:
                st.write("Error:", e)        
    agree_draw_smiles = st.checkbox("Draw chemical structure")
    if agree_draw_smiles:
        smile_code = st_ketcher()
        st.markdown(f"SMILES: {smile_code}")   
    st.write('Serotonergic activity specifically pertains to the functioning of serotonin receptors, integral components of the central nervous system. These receptors play a pivotal role in mediating the effects of serotonin, influencing mood, sleep, and cognitive processes.')
    st.write('Developing a model that discerns the serotoninergic activity of molecules, with a focus on serotonin receptors, is essential for advancing targeted interventions in psychiatric and neurological conditions.')
    st.write('Such a model holds promise in streamlining drug discovery efforts and facilitating personalized approaches by specifically assessing the impact on serotonin receptor function.')
    st.write("Based on the activity model, you can determine whether the molecule is serotonergically active. In the next step, you can check the selectivity to a given serotonin receptor - the **'Selectivity'** tab.")
    st.write("In the case of predictable effective affinity for multiple receptors, one should move on to **'5-HT receptors'** section and, based on these models, check the affinity for all types of serotonin receptors.")
    
    


elif selected == "Selectivity":
    selectivity_model = AutoML(selective_model_path)
    descriptors_for_QSPR = pd.read_table("multiclass_serotonin_receptors_descriptors.txt", header=None)[0][0]
    calc = Calculator(descriptors, ignore_3D=True)
    calc.descriptors = [d for d in calc.descriptors if str(d) in descriptors_for_QSPR] 
    st.markdown('<h1 class="text-second-title">Predictions of selectivity towards serotonin receptors</h1>', unsafe_allow_html=True)
    st.write("Model was built only on active compounds, before using this module check **'Serotonergic activity'**.")
    smiles_input = st.text_input("Input SMILES", key="text")
    col1, col2 = st.columns(2)
    if smiles_input:
        try:
            molecule = Chem.MolFromSmiles(smiles_input)
            if molecule:
                img = Draw.MolToImage(molecule)
                with col1:
                    st.image(img, caption='Chemical structure', use_column_width=True)
            else:
                pass
        except Exception as e:
            st.error(f"Wystąpił błąd: {str(e)}")
    if smiles_input:
        try:
            molecule = Chem.MolFromSmiles(smiles_input)
            if molecule is not None:
                descriptors_value = calc.pandas([molecule])
                descriptors_value_df = pd.DataFrame(descriptors_value)
                for column in descriptors_value_df.select_dtypes(include=['object']).columns:
                    descriptors_value_df[column] = 0
                with col2:
                    with st.spinner('Calculation in progress'):
                        prediction = selectivity_model.predict(descriptors_value_df)
                    def assign_category(value):
                        if value <= 1.5:
                            return '5-HT1A'
                        elif value <= 2.5:
                            return '5-HT1B'
                        elif value <= 3.5:
                            return '5-HT1D'
                        elif value <= 4.5:
                            return '5-HT2A'
                        elif value <= 5.5:
                            return '5-HT2B'    
                        elif value <= 6.5:
                            return '5-HT2C'
                        elif value <= 7.5:
                            return '5-HT3'    
                        elif value <= 8.5:
                            return '5-HT4'
                        elif value <= 9.5:
                            return '5-HT5A' 
                        elif value <= 10.5:
                            return '5-HT6' 
                        else: 
                            return '5-HT7'
                            value <= 11.0
                    prediction_class = assign_category(prediction)
                    st.markdown(f'<div style="{success_style}">Done!</div>', unsafe_allow_html=True)
                    st.write(f"Compound is selective towards **{prediction_class}** serotonin receptor")
                    st.button("Clear SMILES", on_click=clear_text)
                    list_of_important_descriptors = ['SddssS', 'Xch-5d', 'AATSC2s', 'MDEC-33', 'ETA_dPsi_B', 'AATS6s', 'SpMAD_DzZ', 'ATSC3c', 'NsssCH', 
                                 'SaasN']
                    min_values = {'SddssS': -11.181565977471235, 'Xch-5d': 0.0,'AATSC2s': -0.1405448574581022, 'MDEC-33': 0.0004432315575864,
                                  'ETA_dPsi_B': 0.0, 'AATS6s': 1.2037037037037035, 'SpMAD_DzZ': 3.821325332751482,
                                  'ATSC3c': -1.1030119940855148, 'NsssCH': 0, 'SaasN': -0.0249944392971501}
                    max_values = {'SddssS': 0.0, 'Xch-5d': 1.6594659464934112, 'AATSC2s': 2.496406933751847,  'MDEC-33': 27.18089934675436,
                                   'ETA_dPsi_B': 0.0367598784194531,  'AATS6s': 7.930007604222582, 'SpMAD_DzZ': 32.06044269906967,
                                   'ATSC3c': 1.570147655872653, 'NsssCH': 19, 'SaasN': 6.622971712668626}
                    normalized_descriptors_df = (descriptors_value_df - pd.Series(min_values)) / (pd.Series(max_values) - pd.Series(min_values))
                    values_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
                    values_2 = normalized_descriptors_df[list_of_important_descriptors].max().to_list()
                    values_3 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                    labels = normalized_descriptors_df[list_of_important_descriptors].columns
                    desc_condition = sum([val_3 <= val_2 <= val_1 for val_1, val_2, val_3 in zip(values_1, values_2, values_3)])
                    values_1 += values_1[:1]
                    values_2 += values_2[:1]
                    values_2 = [-1 if value < -1 else value for value in values_2]
                    values_2 = [1.5 if value > 1.5 else value for value in values_2]
                    values_3 += values_3[:1]
                    num_labels = len(labels)
                    angles = [n / float(num_labels) * 2 * np.pi for n in range(num_labels)]
                    angles += angles[:1]
                    fig = plt.figure(figsize=(8,8))
                    ax = fig.add_subplot(111, polar=True)
                    color_1 = '#A6A6A6'
                    color_2 = '#4282AA'
                    ax.plot(angles, values_1, color=color_1, label="training set")
                    ax.fill(angles, values_1, alpha=0.25, color=color_1)
                    ax.plot(angles, values_3, color="white")
                    ax.fill(angles, values_3, color='white', alpha=1, edgecolor="white")
                    ax.plot(angles, values_2, color=color_2, label="tested compound", linewidth=3)
                    ax.fill(angles, values_2, alpha=0)
                    ax.set_thetagrids(np.degrees(angles[:-1]), labels)
                    ax.set_ylim(min(min(values_1), min(values_2), min(values_3)), max(max(values_1), max(values_2), max(values_3)))
                    ax.legend()
                    plt.text(0.08, -0.09, "Min-max normalization was applied to descriptors' values based on the training set", ha='left', va='bottom', transform=plt.gca().transAxes, fontsize = 10, color="gray")
                    plt.tight_layout()
                    st.pyplot(fig)    
                    if desc_condition > 6:
                        st.write("Compound is under applicability domain")
                    else:
                        st.write("Compound is not under applicability domain, prediction may be inaccurate")
            else:
                st.write("Invalid SMILES")
        except Exception as e:
            st.write("Error:", e)
    agree_draw_smiles = st.checkbox("Draw chemical structure")
    if agree_draw_smiles:
        smile_code = st_ketcher()
        st.markdown(f"SMILES: {smile_code}")
    st.write('Selectivity towards serotonin receptors is crucial in drug development as it enables the targeted modulation of specific receptor subtypes, minimizing off-target effects and enhancing therapeutic efficacy.')
    st.write('This selectivity allows for a more refined and tailored approach in treating various psychiatric and neurological disorders, such as depression and anxiety, by specifically influencing the serotonin pathways associated with these conditions.')
    st.write('Achieving a high level of selectivity can also reduce the likelihood of adverse reactions, contributing to improved safety profiles and better patient outcomes in pharmacological interventions.')
    st.write("The selectivity model has limitations for compounds that act on more than one serotonin receptor. For a ligand that acts on multiple serotonin receptors, check the affinity value for other types of serotonin receptors - the **'5-HT Receptors'** section.")
             
    
elif selected == "Q&A":
    st.markdown('<h1 class="text-second-title">Questions and answers</h1>', unsafe_allow_html=True)
    st.write('In this part you can find information about possible issues and solutions')
    st.write("---")
    
    with st.container():
        st.subheader(":small_blue_diamond: Calculations do not work or take too long")
        st.write("Sometimes it happens that the application does not work properly. If the calculation is correct, after entering SMILES, a spinning circle will appear with the words 'Calculation in progress.' If this circle does not appear, refresh the application and re-enter the query")
        st.write("---")
    
    with st.container():
        st.subheader(":small_blue_diamond: Invalid SMILES")
        st.write("If the SMILES you enter triggers an 'Invalid SMILES' message, pay attention to whether the entered SMILES is correct. Typically, when the program fails to recognize a molecule, the reason is incorrect capitalization. For example, in the case of carbon, it should be represented as 'C' in uppercase, while chlorine should be represented as 'Cl'.")
        
    with st.container():
        st.write("---")
        st.subheader(":small_blue_diamond: Calculation time")
        st.write("Sometimes the calculation time in single mode or batch mode may be longer than expected. If the calculation for a single molecule takes longer than 1 minute, please refresh the website. For batch mode calculations, the time required depends on the number of SMILES and the number of selected targets, and it can be quite long. The optimal number of molecules for all possible serotonergic targets is 20.")
        
    with st.container():
        st.write("---")
        st.subheader(":small_blue_diamond: Knowledge resources")
        st.write("If you are wondering what is the source of the scientific information presented in this application, it was given based on the following articles:")        
        st.write("- Sharp T, Barnes NM. Central 5-HT receptors and their function; present and future. Neuropharmacology. 2020;177:108155. doi:10.1016/j.neuropharm.2020.108155.")        
        st.write("- Rudnick G, Sandtner W. Serotonin transport in the 21st century. J Gen Physiol. 2019;151(11):1248-1264. doi:10.1085/jgp.201812066.")        
        st.write("- Newman-Tancredi A, Depoortère RY, Kleven MS, Kołaczkowski M, Zimmer L. Translating biased agonists from molecules to medications: Serotonin 5-HT1A receptor functional selectivity for CNS disorders. Pharmacol Ther. 2022;229:107937. doi:10.1016/j.pharmthera.2021.107937.")   
        st.write("- Czub N, Szlęk J, Pacławski A, Klimończyk K, Puccetti M, Mendyk A. Artificial Intelligence-Based Quantitative Structure-Property Relationship Model for Predicting Human Intestinal Absorption of Compounds with Serotonergic Activity. Mol Pharm. 2023;20(5):2545-2555. doi:10.1021/acs.molpharmaceut.2c01117.")
        st.write("- Sivandzade F, Cucullo L. In-vitro blood-brain barrier modeling: A review of modern and fast-advancing technologies. J Cereb Blood Flow Metab. 2018;38(10):1667-1681. doi:10.1177/0271678X18788769.")
        
    with st.container():
        st.write("---")
        st.subheader(":small_blue_diamond: Databases")
        st.write("In this project, we created QSAR and QSPR models based on following databases:")
        st.write("- ZINC https://zinc.docking.org/ ")
        st.write("- ChEMBL https://www.ebi.ac.uk/chembl/")
        st.write("- OCHEM https://ochem.eu/home/show.do")
        st.write("- Data available in this article: Miao, R., Xia, LY., Chen, HH. et al. Improved Classification of Blood-Brain-Barrier Drugs Using Deep Learning. Sci Rep 9, 8802 (2019). https://doi.org/10.1038/s41598-019-44773-4")
    with st.container():
        st.write("---")
        st.subheader(":small_blue_diamond: Models' structure")
        st.write("For QSAR and QSPR models creation we used **mljar-supervised** tool (https://supervised.mljar.com/).")
        st.write("Here are presented detailed structure of each model:")
        st.write("- **5-HT1A** - LightGBM (last version Sept. 2023)")
        st.write("- **5-HT1B** - Ensemble model (5 x Xgboost, 4 x LightGBM) (last version Sept. 2023)")
        st.write("- **5-HT1D** - Ensemble model (5 x Xgboost, 4 x LightGBM, 2 x CatBoost) (last version Sept. 2023)")
        st.write("- **5-HT2A** - Xgboost (last version Sept. 2023)")
        st.write("- **5-HT2B** - Ensemble model (6 x Xgboost, 2 x LightGBM, 5 x CatBoost) (last version Sept. 2023)")
        st.write("- **5-HT2C** - LightGBM (last version Sept. 2023)")
        st.write("- **5-HT3** - Ensemble model (2 x Xgboost, 3 x CatBoost, 1 x LightGBM) (last version Sept. 2023)")
        st.write("- **5-HT4** - Xgboost (last version Nov. 2023)")
        st.write("- **5-HT5A** - Ensemble model (3 x Xgboost, 3 x CatBoost, 1 x EstraTrees, 1 x LightGBM) (last version Sept. 2023)")
        st.write("- **5-HT6** - Xgboost (last version Sept. 2023)")
        st.write("- **5-HT7** - CatBoost (last version Sept. 2023)")
        st.write("- **SERT** - LightGBM (last version Sept. 2023)")
        st.write("- **BBB** - Ensemble model (5 x Xgboost, 7 x NeuralNetwork, 2 x LightGBM) (last version Sept. 2023)")
        st.write("- **HIA** - AI-based system [Classification: Ensemble model (3 x CatBoost, 1 x RandomForest, 1 x ExtraTrees, 2 x NeuralNetwork); Regression: Ensemble model (5 x RandomForest, 2 x NeuralNetwork)]. AI-based system have been already published: Czub N, Szlęk J, Pacławski A, Klimończyk K, Puccetti M, Mendyk A. Artificial Intelligence-Based Quantitative Structure-Property Relationship Model for Predicting Human Intestinal Absorption of Compounds with Serotonergic Activity. Mol Pharm. 2023;20(5):2545-2555. https://pubs.acs.org/doi/10.1021/acs.molpharmaceut.2c01117.")
        st.write("- **Serotonergic activity** - Ensemble model (5 x Xgboost, 6 x LightGBM, 2 x NeuralNetwork, 2 x CatBoost) (last version Jan. 2024)")
        st.write("- **Selectivity** - Ensemble model (1 x Baseline, 5 x LightGBM, 5 x DecisionTree, 2 x Xgboost) (last version Jan. 2024)")
    with st.container():
        st.write("---")
        st.subheader(":small_blue_diamond: SMILES type")
        st.write("During SerotoninAI usage enter SMILES without making any changes, as the models have been trained based on the basic SMILES representation found in databases such as DrugBank, ChEMBL, and ZINC.")

elif selected == "Contact":
    st.markdown('<h1 class="text-second-title">Contact</h1>', unsafe_allow_html=True)
    st.write("My name is **Natalia Łapińska** (maiden name Czub) and I am the author of ***SerotoninAI***. I hope this application will help develop therapies for new drugs affecting the central nervous system. ")
    st.write('A few words about me, I am a Polish woman 🇵🇱 who fell in love with the world of artificial intelligence and the positive opportunities it gives us to improve all areas of research.')
    st.write("Wanting to define myself, one term is not enough. Here's a list of the ones that best fit scientist :microscope:, pharmacist 💊, almost PhD 👩‍🎓, Python programmer :snake:, cheminformatics :female-technologist:, data scientist 📊 and AI specialist🔥.")
    st.write("If you want to develop AI-based drug discovery together with me, I invite you to contact me!")
    st.markdown(""" <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.3/css/all.min.css" rel="stylesheet"> """, unsafe_allow_html=True)
    st.markdown('<i class="fab fa-linkedin"></i>  www.linkedin.com/in/natalia-czub', unsafe_allow_html=True)
    st.markdown('<i class="fab fa-github-square"></i>  https://github.com/nczub', unsafe_allow_html=True)
    st.markdown('<i class="fa-solid fa-message"></i>  lapinska.natalia@inbox.eu', unsafe_allow_html=True)
