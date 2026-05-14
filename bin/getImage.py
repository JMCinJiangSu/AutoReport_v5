#-*- coding:gbk -*-
from docxtpl import InlineImage
from docx.shared import Mm
import os
import re
from libs.getImageSize import image_file

'''
Discription
	
	该脚本用来处理填充的图片，主要处理使用插入方法（InlineImage）的图片

'''
def render_image(tpl, data, jsonDict, report_name, image, config):
	result = {}
	# 插入图片
	# 1. BRCA IGV图
	result["brca_all"] = [InlineImage(tpl, i, width=Mm(150)) for i in data["var_brca"]["igvplot"]["all"] if i != None and os.path.exists(i)]
#	print (result["brca_all"])
	result["brca_g"] = [InlineImage(tpl, i, width=Mm(150)) for i in data["var_brca"]["igvplot"]["germline"] if i != None and os.path.exists(i)]
	# 2. MSI图
	if data["msi"] and "img_path" in data["msi"].keys() and data["msi"]["img_path"] and os.path.exists(data["msi"]["img_path"]):
		result["msi"] = InlineImage(tpl, data["msi"]["img_path"], width=Mm(90))
		# MSI加个尺寸小点的图-2024.09.26
		result["msi_bds"] = InlineImage(tpl, data["msi"]["img_path"], width=Mm(70))
	# 3. TME图
	#    新加个小点的尺寸-适用江门中心MP-2024.04.22
	if data["qc"] and "rna_data_qc" in data["qc"].keys() and data["qc"]["rna_data_qc"] and "tmeplot" in data["qc"]["rna_data_qc"].keys() and \
		data["qc"]["rna_data_qc"]["tmeplot"] and os.path.exists(data["qc"]["rna_data_qc"]["tmeplot"]):
		result["tme"] = InlineImage(tpl, data["qc"]["rna_data_qc"]["tmeplot"], width=Mm(110))
		result["tme_jmxz"] = InlineImage(tpl, data["qc"]["rna_data_qc"]["tmeplot"], width=Mm(90))
	# 2025.06.26-TME图片存放路径改了
	elif "tme_type" in jsonDict.keys() and jsonDict["tme_type"] and "tmeplot" in jsonDict["tme_type"].keys() and jsonDict["tme_type"]["tmeplot"]:
		# 系统返回json这个路径不是完整的，先拼接下，等报告系统加好再改
		# v3
		#tmeplot_path = os.path.join("/var/www/html/report_backend/system_files/json_file/", jsonDict["tme_type"]["json_batch_name"], jsonDict["tme_type"]["tmeplot"])
		# v4
		tmeplot_path = os.path.join("/var/www/html/report_backend/system_files/box/tmp/un7z", jsonDict["tme_type"]["json_batch_name"], jsonDict["tme_type"]["tmeplot"])
		if os.path.exists(tmeplot_path):
			result["tme"] = InlineImage(tpl, tmeplot_path, width=Mm(110))
			result["tme_jmxz"] = InlineImage(tpl, tmeplot_path, width=Mm(90))
	# 2025.06.26-更新完成
	# 4. GEP图 先判断gep.img_path，有的话直接展示，无则判断qc.rna_data_qc.gepplot
	#    新加个小点的尺寸-适用江门中心MP-2024.04.22
	if data["gep"] and "img_path" in data["gep"].keys() and data["gep"]["img_path"] and os.path.exists(data["gep"]["img_path"]):
		result["gep"] = InlineImage(tpl, data["gep"]["img_path"], width=Mm(90))
		result["gep_jmzx"] = InlineImage(tpl, data["gep"]["img_path"], width=Mm(75))
	elif data["qc"] and "rna_data_qc" in data["qc"].keys() and data["qc"]["rna_data_qc"] and "gepplot" in data["qc"]["rna_data_qc"].keys() and \
		data["qc"]["rna_data_qc"]["gepplot"] and os.path.exists(data["qc"]["rna_data_qc"]["gepplot"]):
		result["gep"] = InlineImage(tpl, data["qc"]["rna_data_qc"]["gepplot"], width=Mm(90))
		result["gep_jmzx"] = InlineImage(tpl, data["qc"]["rna_data_qc"]["gepplot"], width=Mm(75))
	# 2025.08.20-GEP新增路径
	elif data["gep"] and "json_batch_name" in data["gep"].keys() and data["gep"]["json_batch_name"] and "gepplot" in data["gep"].keys() and data["gep"]["gepplot"]:
		gepplot_path = os.path.join("/var/www/html/report_backend/system_files/box/tmp/un7z", data["gep"]["json_batch_name"], data["gep"]["gepplot"])
		if os.path.exists(gepplot_path):
			result["gep"] = InlineImage(tpl, gepplot_path, width=Mm(90))
			result["gep_jmzx"] = InlineImage(tpl, gepplot_path, width=Mm(75))
	# 2025.08.20-新增完成
	# 5. 新增I/II类变异Snvindel IGV图-2022.10.20
	result["igv_I_II"] = [InlineImage(tpl, i, width=Mm(180)) for i in data["var"]["igv_I_II"] if os.path.exists(i)]
	# 新增适用福建附一CP40的IGV图，长11cm，宽5cm-2023.10.20
	# 加一个限制，IGV图格式不能是svg的-2024.08.14
	result["igv_I_II_FJFY"] = [InlineImage(tpl, i, width=Mm(110), height=Mm(50)) for i in data["var"]["igv_I_II"] if re.split("\.", i)[-1] != "svg" and os.path.exists(i)]
	# 6. 新增4/5类变异Snvindel IGV图-2022.11.18
	result["igv_4_5"] = [InlineImage(tpl, i, width=Mm(180)) for i in data["var"]["igv_4_5"] if os.path.exists(i)]
	# 新增适用福建附一HRR的IGV图，长11cm，宽5cm-2024.05.15
	result["igv_4_5_FJFY"] = [InlineImage(tpl, i, width=Mm(110), height=Mm(50)) for i in data["var"]["igv_4_5"] if os.path.exists(i)]
	# 7. TMB图
	if data["tmb"] and "img_path" in data["tmb"].keys() and os.path.exists(data["tmb"]["img_path"]):
		result["tmb"] = InlineImage(tpl, data["tmb"]["img_path"], width=Mm(100))
		# TMB加个尺寸小点的图-2024.09.26
		result["tmb_bds"] = InlineImage(tpl, data["tmb"]["img_path"], width=Mm(80))
	# 8. PD-L1图
	if data["pdl1"] and "file_pdl1" in data["pdl1"].keys() and os.path.exists(data["pdl1"]["file_pdl1"]):
		result["pdl1"] = InlineImage(tpl, data["pdl1"]["file_pdl1"], width=Mm(80))
	# 9. CNV图
	if "cnv_file_path" in jsonDict.keys() and jsonDict["cnv_file_path"] and "abs_path" in jsonDict["cnv_file_path"].keys() and \
		jsonDict["cnv_file_path"]["abs_path"] and os.path.exists(jsonDict["cnv_file_path"]["abs_path"]):
		result["cnv"] = InlineImage(tpl, jsonDict["cnv_file_path"]["abs_path"], width=Mm(160))
	# 10. MLPA图（仅del）
	result["mlpa_image_del"] = [InlineImage(tpl, i, width=Mm(150)) for i in data["var_brca"]["mlpa_image_del"] if os.path.exists(i)]

	# 新增肿瘤发生发展相关Snindel IGV图-2023.06.27
	result["igv_onconodrug"] = [InlineImage(tpl, i, width=Mm(180)) for i in data["var"]["igv_onconodrug"] if os.path.exists(i)]
	# 新增结束-2023.06.27

	# 新增福建附一肿瘤发生发展相关/III类变异IGV图-2025.08.08
	result["igv_III_onconodrug_FJFY"] = [InlineImage(tpl, i, width=Mm(110), height=Mm(50)) for i in data["var"]["igv_III_onconodrug_FJFY"] if os.path.exists(i)]
	# 2025.08.08-添加完成

	# 新增福建附一150 林奇5基因变异IGV图-2025.08.08
	result["igv_4_5_lyn"] = [InlineImage(tpl, i, width=Mm(110), height=Mm(50)) for i in data["var"]["igv_4_5_lyn"] if os.path.exists(i)]
	# 2025.08.08-添加完成

	# 获取配置表中图片指定尺寸，默认为100
	image_sise_dict = image_file(config)
	# 模板中固定的图片，直接拼接会显示不出来，这边把固定图片按插入的形式进行填充
	image_list = [i for i in os.listdir(image)]
	for i in image_list:
		image_name = re.split("\.", i)[0]
		width = image_sise_dict.get(image_name, 100) 
		result["fixed_"+image_name] = InlineImage(tpl, image+"/"+i, width=Mm(width))
	
	# 新增cnv图-适用CP40-2024.08.07
	# cp200 cnv图从qc移动到cnv_file_path，这边做兼容-2025.04.18
	cnvplot_path = ""
	if "dna_data_qc" in data["qc"].keys() and data["qc"]["dna_data_qc"] and "cnvplot" in data["qc"]["dna_data_qc"].keys() and data["qc"]["dna_data_qc"]["cnvplot"]:
		# 系统返回json这个路径不是完整的，先拼接下，等报告系统加好再改
		# v3
		#cnvplot_path = os.path.join("/var/www/html/report_backend/system_files/json_file/", data["sample"]["json_batch_name"], data["qc"]["dna_data_qc"]["cnvplot"])
		# v4
		cnvplot_path = os.path.join("/var/www/html/report_backend/system_files/box/tmp/un7z", data["sample"]["json_batch_name"], data["qc"]["dna_data_qc"]["cnvplot"])
	elif "cnv_file_path" in jsonDict.keys() and jsonDict["cnv_file_path"] and "abs_path" in jsonDict["cnv_file_path"].keys() and jsonDict["cnv_file_path"]["abs_path"]:
		cnvplot_path = jsonDict["cnv_file_path"]["abs_path"]
	# 阴性的没有返回cnv图路径，这里拼接，过渡使用-2025.04.25
	# sample_id = sample_id+library_id+flowcell_lane
	# /var/www/html/report_backend/system_files/json_file/json_batch_name/sample_id/sample_id.cnv.png 
	else:
		if jsonDict["qc"] and "dna_data_qc" in jsonDict["qc"].keys() and jsonDict["qc"]["dna_data_qc"] and type(jsonDict["qc"]["dna_data_qc"]).__name__=="list":
			sample_id = jsonDict["sample_info"]["sample_id"]
			library_id = jsonDict["qc"]["dna_data_qc"][0]["library_id"] if "library_id" in jsonDict["qc"]["dna_data_qc"][0].keys() and jsonDict["qc"]["dna_data_qc"][0]["library_id"] else ""
			flowcell_lane = jsonDict["qc"]["dna_data_qc"][0]["flowcell_lane"] if "flowcell_lane" in jsonDict["qc"]["dna_data_qc"][0].keys() and jsonDict["qc"]["dna_data_qc"][0]["flowcell_lane"] else ""
			cnv_sample_id = "{0}_{1}_{2}".format(sample_id, library_id, flowcell_lane)
			# v3
			#cnvplot_path = os.path.join("/var/www/html/report_backend/system_files/json_file/", data["sample"]["json_batch_name"], cnv_sample_id, cnv_sample_id+".cnv.png")
			# v4
			cnvplot_path = os.path.join("/var/www/html/report_backend/system_files/box/tmp/un7z", data["sample"]["json_batch_name"], cnv_sample_id, cnv_sample_id+".cnv.png")
	print (cnvplot_path)
	# 2025.04.25-更新完成

		#cnvplot_path = data["qc"]["dna_data_qc"]["cnvplot"]
	if os.path.exists(cnvplot_path):
		result["cnvplot"] = InlineImage(tpl, cnvplot_path, width=Mm(110), height=Mm(80.5))
		# 2025.02.14-新增一个福建附一CP200图，要求尺寸高25.35cm，宽18.12cm
		result["cnvplot_cp200_FJFY"] = InlineImage(tpl, cnvplot_path, width=Mm(181.2), height=Mm(253.5))
		# 2025.02.14-新增完成
		# 2025.04.10-新增福建省立CP200图，要求尺寸高21.7cm，宽17.75cm
		result["cnvplot_cp200_FJSL"] = InlineImage(tpl, cnvplot_path, width=Mm(177.5), height=Mm(217))
		# 2025.04.10-新增完成
		# 2025.04.25-cp200新流程cnv图大小有改动-新增一个字段
		result["cnvplot_v2"] = InlineImage(tpl, cnvplot_path, width=Mm(180.1), height=Mm(117.4))
		# 2025.04.25-新增完成
	# 新增结束-2024.08.07

	# 新增华西送检信息照片和条码号-2025.05.19
	# 2025.07.08-华西照片和条码号路径需要拼接
	hx_v4_path = "/var/www/html/report_backend/system_files/images"
	patient_photo_path = jsonDict["sample_info"]["patient_photo"] if "patient_photo" in jsonDict["sample_info"].keys() and jsonDict["sample_info"]["patient_photo"] else ""
	hx_patient_photo_path = os.path.join(hx_v4_path, patient_photo_path)
	#if os.path.exists(patient_photo_path):
	if os.path.exists(hx_patient_photo_path) and patient_photo_path:
		result["patient_photo"] = InlineImage(tpl, patient_photo_path, height=Mm(35))
	barcode_photo_path = jsonDict["sample_info"]["barcode_photo"] if "barcode_photo" in jsonDict["sample_info"].keys() and jsonDict["sample_info"]["barcode_photo"] else ""
	hx_barcode_photo_path = os.path.join(hx_v4_path, barcode_photo_path)
	#if os.path.exists(barcode_photo_path):
	if os.path.exists(hx_barcode_photo_path) and barcode_photo_path:
		result["barcode_photo"] = InlineImage(tpl, barcode_photo_path, width=Mm(50))
	# 2025.05.19-新增完成
	# 2025.05.22-新增华西，不存在照片和条码图片的话，插入下面的预留字段
	result["barcode_photo_raw"] = "{{barcode_photo}}"
	result["patient_photo_raw"] = "{{patient_photo}}"
	# 2025.05.22-新增完成

	# 2025.12.02-新增MRD图片
	if data["mrd"]["mrd_last"] and "mrd_img_path" in data["mrd"]["mrd_last"].keys() and data["mrd"]["mrd_last"]["mrd_img_path"] and \
		os.path.exists(data["mrd"]["mrd_last"]["mrd_img_path"]):
		result["mrd"] = InlineImage(tpl, data["mrd"]["mrd_last"]["mrd_img_path"], height=Mm(145))
	# 2025.12.02-新增完成

	# 插入或替换图片
	# 7. TMB图
#	if data["tmb"] and "img_path" in data["tmb"].keys() and os.path.exists(data["tmb"]["img_path"]):
#		# 报告打开报错，TMB图插入方式修改-仅临检通用Master样本（2022.09.27新增复旦中山）
#		if re.search("Master|90|87", data["sample"]["prod_names"]) and re.search("rummage|FDZS|SDSL", report_name):
#			tpl.replace_pic("test_TMB.jpg", data["tmb"]["img_path"])
#			print ("replace")
#		else:
#			result["tmb"] = InlineImage(tpl, data["tmb"]["img_path"], width=Mm(100))
#			print (inline)

	# 替换图片
	# 8. PD-L1图 图片太多可能会报错，PDL1用替换
#	if data["pdl1"] and "file_pdl1" in data["pdl1"].keys() and os.path.exists(data["pdl1"]["file_pdl1"]):
#		tpl.replace_pic("test_PDL1.jpg", data["pdl1"]["file_pdl1"])
	# 9. CNV图
	# 曲线贴图-这个原因待查找！！！
	# 使用插入图片的方法，本地测试正常、线上运行打开word会报错，可能的原因有1）word格式问题，2）模块版本差异；排查太耗时了，暂时放弃
	# 改用replace方法，插入图片名跟本地不符合，代码找不到，排查耗时，先不找了。暂时使用之前正常替换图片的图片模板，成功运行，先用着吧，有空再来排查
#	if re.search("Pan116（组织）|LC76（组织）|CRC25（组织）|GA18（组织）|TC21（组织）", data["sample"]["prod_names"]) and re.search("rummage", report_name):
#		if "cnv_file_path" in jsonDict.keys() and jsonDict["cnv_file_path"] and "abs_path" in jsonDict["cnv_file_path"].keys() and jsonDict["cnv_file_path"]["abs_path"] and os.path.exists(jsonDict["cnv_file_path"]["abs_path"]):
#			tpl.replace_pic("test_MSI.png", jsonDict["cnv_file_path"]["abs_path"])
	
	return result
