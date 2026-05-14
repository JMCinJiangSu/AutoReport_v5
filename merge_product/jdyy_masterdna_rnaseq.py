#-*- coding:gbk -*-
import os
import json
from bin import getSampleInfo
from bin import getQC
from bin import getVar
from bin import getRNAexp
from bin import getChemo
from bin import getGEP
from bin import getTME
from bin import getRefence
from bin import getApprovalRegimen
from bin import getClinicTrial
from bin import getMSI
from bin import getTMB
from bin import getHRD
from bin import getPDL1
from bin import getMRD

#---------------------
# 本程序用来处理吉大一MasterDNA+RNAseq出一份报告的需求
# 由于MP没有RNA的结果，可将RNAseq的结果直接传入MP中，然后按单个样本进行数据转化
#----------------------

config = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "config")

# 1. 遍历样本信息，逐个样本进行处理，将RNAseq的原始数据整合到MPDNA里
def stran_json_data(relation_jsonDict, jsonDict_input):
    MP_json = {}
    RNA_json = {}
    if jsonDict_input["sample_info"]["raw_prod_names"] == "Master Panel（组织）":
        MP_json = jsonDict_input
    else:
        MP_json = relation_jsonDict
    
    if jsonDict_input["sample_info"]["raw_prod_names"] == "RNASeq（肉瘤）":
        RNA_json = jsonDict_input
    else:
        RNA_json = relation_jsonDict
    MP_json["sample_info"]["rna_sample_id"] = RNA_json["sample_info"]["rna_sample_id"]
    MP_json["qc"]["rna_data_qc"] = RNA_json["qc"]["rna_data_qc"]
    MP_json["lib_quality_control"]["rna_lib_qc"] = RNA_json["lib_quality_control"]["rna_lib_qc"]
    MP_json["rna_sv"] = RNA_json["rna_sv"]
    MP_json["refer"].extend(RNA_json["refer"])

    return MP_json

# 2. 主函数-处理后的数据进行转化，返回合并后的数据
def MergeReport_main(relation_jsonDict, jsonDict_input, report_name):
    data = {}
    if relation_jsonDict and jsonDict_input:
        jsonDict = stran_json_data(relation_jsonDict, jsonDict_input)
        jsonDict["sample_info"]["tumor_names_cn"] = jsonDict["sample_info"]["tumor_type"]
        data["sample"] = getSampleInfo.getSample(jsonDict)
        data["qc"], data["lib_quality_control"] = getQC.getJsonQC(jsonDict)
        data["therapeutic_regimen"] = getApprovalRegimen.getRegimen(jsonDict)
        data["clinic_trial"] = getClinicTrial.getClinic(jsonDict)
        data["msi"] = getMSI.getMSI(jsonDict, config)
        data["tmb"] = getTMB.getTMB(jsonDict, config)
        data["var"], data["var_brca"], data["var_hrr_shsy"] = getVar.getVar(jsonDict, config, report_name)
        data["rna_exp"] = getRNAexp.getRNA_exp(jsonDict)
        data["chemo"] = getChemo.getchemo(jsonDict, config)
        data["gep"] = getGEP.getGEP(jsonDict)
        data["tme"], data["tme_score"] = getTME.getTME(jsonDict)
        data["hrd"] = getHRD.getHRD(jsonDict, data["var"]["ec_type"]["BRCA1_level12"]+data["var"]["ec_type"]["BRCA2_level12"], config)
        data["refer"] = {}
        data["refer"]["fixed"] = getRefence.getfixed_refer(report_name, data["sample"]["tumor_list"], data["var"]["knb"], data["msi"], config)
        data["refer"]["dynamic"] = getRefence.getdynamic_refer(jsonDict, data["var"], data["msi"], data["hrd"], data["var_brca"])
        data["pdl1"] = getPDL1.getPDL1(jsonDict)
        data["mrd"] = getMRD.getMRD(jsonDict)

    return data