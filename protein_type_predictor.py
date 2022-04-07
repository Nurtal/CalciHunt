


def get_target_list():
    """
    return the list of all possible protein type 1
    found in dataset
    """

    ## importation
    import pandas as pd

    ## parameters
    target_list = []

    ## load data
    df = pd.read_csv("calcibot_dataset.csv")
    for target in df["protein type 1"]:
        if(target not in target_list):
            target_list.append(target)
            print(target)

    ## return target list
    return target_list



def get_summary(gene):
    """
    """

    ## importation
    import mygene
    import pprint

    ## init session
    mg = mygene.MyGeneInfo()

    ## get gene id
    stuff = mg.query('symbol:'+str(gene), species='human')
    gene_id = stuff["hits"][0]["_id"]

    ## extract summary
    summary =  mg.getgene(gene_id, fields='summary')
    summary = summary['summary']

    ## return summary
    return summary


def craft_train_dataset():
    """
    """

    ## importation
    import pandas as pd

    ## parameters
    output_file_name = "protein_type_1_train_dataset.csv"

    ## load dataset
    df = pd.read_csv("calcibot_dataset.csv")

    ## init output
    output_data = open(output_file_name, "w")
    output_data.write("GENE\tSUMMARY\tCLASS\n")

    ## loop over row
    for index, row in df.iterrows():

        gene = row['Gene name']
        protein = row["protein type 1"]
        try:
            summary = get_summary(gene)
            summary = summary.replace("\n", "")
            summary = summary.replace("\t", "")
        except:
            summary = "nan"

        ## update file
        output_data.write(str(gene)+"\t"+str(summary)+"\t"+str(protein)+"\n")

        ## display something
        print("[+] -> "+str(gene))
        print("[+] "+str(summary))
        print("----------------")

    ## close file
    output_data.close()


def hunt_cytoskeleton_prot(summary):
    """
    """

    ## importation
    import re

    ## parameters
    target_list = ["cytoskeleton", "cytoskeletal", "structural component", "actin filaments"]
    hit = False

    ## hunt target
    for target in target_list:
        if(re.search(target, summary)):
            hit = True

    ## return results
    return hit


def evaluate():
    """
    """

    ## importation
    import pandas as pd

    ## parameters
    skeleton_match = 0
    skeleton_miss = 0
    skeleton_falsepos = 0
    skeleton_trueneg = 0

    ## load data
    df = pd.read_csv("protein_type_1_train_dataset.csv", sep="\t", encoding = 'cp1252')
    df = df.dropna()

    ## loop over data
    for index, row in df.iterrows():

        ## extract data
        summary = row["SUMMARY"]
        label = row["CLASS"]

        ## deal with cytoskeleton
        predict_as_skeleton = hunt_cytoskeleton_prot(summary)
        if(label == "Cytoskeleton protein" and predict_as_skeleton):
            skeleton_match +=1
        elif(label =="Cytoskeleton protein" and not predict_as_skeleton):
            skeleton_miss +=1
        elif(label !="Cytoskeleton protein" and predict_as_skeleton):
            skeleton_falsepos +=1
        elif(label !="Cytoskeleton protein" and not predict_as_skeleton):
            skeleton_trueneg +=1
            print(summary)



    ## global eavaluation
    #-> skeleton
    skeleton_acc = (float(skeleton_match)+float(skeleton_trueneg)) / df.shape[0]
    print(skeleton_acc*100.0)



evaluate()
