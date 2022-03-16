




def gene_card_request(ensgid):
    '''
    Passing Ensembl gene id, return its webpage on Genecards
    '''

    ## importation

    import requests
    import html
    import time
    import argparse

    url_to_request='https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + ensgid
    headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36'}
    r = requests.get(url_to_request, headers=headers)
    gene_card_html = html.unescape(r.text)

    return gene_card_html


"""
stuff = gene_card_request("GC17P076453")
print(stuff)
"""

def hunt(genecard_id):
    """
    """

    ## importation
    import re
    from bs4 import BeautifulSoup

    ## parameters
    best_confidence = 0
    best_spot = "NA"

    ## get html page for gene
    html_response = gene_card_request(genecard_id)

    ## hunt location in html
    soup = BeautifulSoup(html_response, 'html.parser')
    stuff = soup.find("div", {"id": "jensenLocalization"})

    extracted = re.findall("(<div alt=\"Jensen Localization Image.+ title=\"Jensen Localization Image)", str(stuff))
    if(len(extracted) > 0):
        extracted = extracted[0]

        ## clean
        extracted = extracted.split("data-model=\"")
        extracted = extracted[1]
        extracted = extracted.replace("{", "")
        extracted = extracted.replace("[", "")
        extracted = extracted.replace("\"", "")
        extracted = extracted.replace(" ", "")
        extracted = extracted.replace("gc-compartment-image=", "")
        extracted = extracted.replace("id=jensenLocalizationImage", "")
        extracted = extracted.replace("src=/Images/jensen_images.html", "")
        extracted = extracted.replace("title=Jensen Localization Image", "")
        extracted = extracted.replace("name:golgi=title=JensenLocalizationImage", "")
        extracted = extracted.replace("}]=title=JensenLocalizationImage", "")
        extracted = extracted.split("},")

        ## extract informations
        spot_to_confidence = {}
        for elt in extracted:
            elt = elt.split(",")
            if(len(elt) > 1):
                spot = elt[0]
                spot = spot.replace("name:", "")
                spot = spot.split("=")
                if(len(spot)>1):
                    spot = spot[1]
                else:
                    spot = spot[0]
                confidence = elt[1]
                confidence = confidence.replace("confidence:", "")

                #-> clean confidence
                if("}" in confidence):
                    confidence = confidence.split("}")
                    confidence = confidence[0]

                #-> update data structure
                spot_to_confidence[spot] = confidence

                #-> pick best spot
                if(float(confidence) > best_confidence):
                    best_confidence = float(confidence)
                    best_spot = spot
                elif(float(confidence) == best_confidence):
                    best_spot = best_spot+" / "+spot

    ## return best confidence
    return best_spot


def run_hunt_on_data_file():
    """
    """
    ## importation
    import pandas as pd
    import os

    ## parameters
    target_file = "localisation_dataset.csv"
    output_file = "localization_hunt_results.csv"
    already_computed = []

    ## stop and go
    if(os.path.isfile(output_file)):
        df = pd.read_csv(output_file, encoding = 'cp1252')
        already_computed = list(df['ID'])
        output_data = open(output_file, "a")
    else:
        #-> init output file
        output_data = open(output_file, "w")
        output_data.write("ID,LOC1,LOC2,LOC3,PREDICTED\n")

    ## load dataset
    df = pd.read_csv(target_file)
    nb_to_process = int(df.shape[0])

    ## loop over data
    cmpt = 0
    for index, row in df.iterrows():

        gen_id = row['Gene name']
        true_spot_1 = row['localization 1']
        true_spot_2 = row['localization 2']
        true_spot_3 = row['localization 3']

        ## hunt spot
        if(gen_id not in already_computed):

            ## hunt
            predicted_spot = hunt(gen_id)

            ## compute progress
            cmpt +=1
            progress = (float(cmpt) / float(nb_to_process))*100.0
            progress += "%"

            ## display something
            print("[+][HUNT] "+str(gen_id)+" => "+str(predicted_spot)+" # "+str(progress)+" completed")

            ## update result file
            output_data.write(str(gen_id)+","+str(true_spot_1)+","+str(true_spot_2)+","+str(true_spot_3)+","+str(predicted_spot)+"\n")
        else:

            ## compute progress
            cmpt +=1
            progress = (float(cmpt) / float(nb_to_process))*100.0
            progress += "%"

            ## display something
            print("[+][HUNT] "+str(gen_id)+" => already extracted"+" # "+str(progress)+" completed")


    ## close output_file
    output_data.close()

run_hunt_on_data_file()
