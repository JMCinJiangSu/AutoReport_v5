#-*- coding:gbk -*-

# 2026.01.06-临检捕获样本加个规则报告中填充的文库总量字段是“library_qty”，需要改为预文库总量
# 临时使用代码替换，后续再慢慢修改模板
# ***文库总量：需要把dna_pre_library_qty的值赋给library_qty
# ***DNA预文库：需要把dna_pre_library_qty的值赋给dna.library_qty
# DNA终文库：dna_fnl_library_qty（oncopro血液使用，无需修改）
# ***RNA预文库：rna_pre_library_qty（需要把rna_pre_library_qty的值赋给rna.library_qty）
# RNA终文库：rna_fnl_library_qty（暂未使用，无需修改）
# 对照样本文库总量：control_library_qty（暂未使用，无需修改）
# handle项目，文库总量使用library_qty（无需修改）
def clinical_libraty_qty(report_module_type, prod_names, lib_data):
    capture_prod = ["10基因（血液）", "10基因（组织）", "CRC25（血液）", "CRC25（组织）", "GA18（血液）", "GA18（组织）", \
                    "LC76（血液）", "LC76（组织）", "Master Panel（血液）", "Master Panel（血液-综合）", "Master Panel（组织）", \
                    "Master Panel（组织_北医）", "Master Panel（组织-综合）", "MRD（LC10）", "OncoPro Panel（血液）-5基因", "OncoPro（血液）", \
                    "Pan116（血液）", "Pan116（血液）-黑色素瘤", "Pan116（血液）-乳腺", "Pan116（血液）-头颈", "Pan116（血液）-胰腺", \
                    "Pan116（组织）", "Pan116（组织）-黑色素瘤", "Pan116（组织）-乳腺", "Pan116（组织）-头颈", "Pan116（组织）-胰腺", \
                    "RNASeq（肉瘤）", "TC21（血液）", "TC21（组织）"]
    if report_module_type == "rummage" and prod_names in capture_prod:
        if lib_data and "lib_dna_qc" in lib_data.keys() and lib_data["lib_dna_qc"] and "dna_pre_library_qty" in lib_data["lib_dna_qc"].keys() and lib_data["lib_dna_qc"]["dna_pre_library_qty"]:
            lib_data["lib_dna_qc"]["library_qty"] = lib_data["lib_dna_qc"]["dna_pre_library_qty"]
            lib_data["lib_dna_qc"]["library_qty_num"] = lib_data["lib_dna_qc"]["dna_pre_library_qty_num"]
            print ("临检", prod_names, "dna_pre_library_qty值", lib_data["lib_dna_qc"]["dna_pre_library_qty"])
        if lib_data and "rna_lib_qc" in lib_data.keys() and lib_data["rna_lib_qc"] and "rna_pre_library_qty" in lib_data["rna_lib_qc"].keys() and lib_data["rna_lib_qc"]["rna_pre_library_qty"]:
            lib_data["rna_lib_qc"]["library_qty"] = lib_data["rna_lib_qc"]["rna_pre_library_qty"]
            lib_data["rna_lib_qc"]["library_qty_num"] = lib_data["rna_lib_qc"]["rna_pre_library_qty_num"]
            print ("临检", prod_names, "rna_pre_library_qty值", lib_data["rna_lib_qc"]["rna_pre_library_qty"])

    return lib_data