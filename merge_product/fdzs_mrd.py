
import re
from bin import getMRD, getQC, getRefence, getSampleInfo
from libs.rule import decimal_percen
from docx import Document
from docx.oxml import OxmlElement
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.shared import Pt

def get_data(jsonDict, outfile, report_name, config):
    data = {}
    data["mrd"] = getMRD.getMRD(jsonDict, outfile)
    data["sample"] = getSampleInfo.getSample(jsonDict)
    data["qc"] = getQC.getJsonQC(jsonDict)[0]
    data["lib_quality_control"] = getQC.getJsonQC(jsonDict)[1]
    data["refer"] = {}
    data["refer"]["fixed"] = getRefence.getfixed_refer(report_name, data["sample"]["tumor_list"], data["var"]["knb"], data["msi"], config)
    data["refer"]["dynamic"] = getRefence.getdynamic_refer(jsonDict, data["var"], data["msi"], data["hrd"], data["var_brca"])
    data["table_data"] = process_merge_order(jsonDict)

# 处理报告json中merge_order内容，主要获取基因、变异、频率、时间
# 针对MRD WES，只需要SNV和SV的内容

def process_merge_order(jsonDict):
    merge_order = jsonDict["merge_order"]
    if not merge_order:
        return []
    snv = merge_order.get("snv")["table_data"] if merge_order.get("snv") else []
    sv = merge_order.get("sv")["table_data"] if merge_order.get("sv") else []

    sample_id_list = merge_order.get("snv")["sample_id_list"] if merge_order.get("snv") else []
    sample_id_list = [item["sample_id"] for item in sample_id_list] if sample_id_list else []
    # "sample_id_list" 缺一个时间戳，按这个时间先后排序，用于确定频率输出先后顺序
    # 区分样本类型
    tissue_sample = [item for item in sample_id_list if item.get("sample_type") == "Tissue"]
    blood_sample = [item for item in sample_id_list if item.get("sample_type") == "Plasma"]

    if len(tissue_sample) > 1:
        raise ValueError("样本中存在多个组织样本，请检查merge_order的sample_id_list内容")
    tissue_sample = tissue_sample[0] if tissue_sample else ""

    headers = ["基因", "检测变异", "组织丰度"]
    for sample in blood_sample:
        sample_id = sample["sample_id"]
        receive_date = sample["receive_date"].split()[0].replace("-", ".") if sample.get("receive_date") else ""
        headers.append(f"血检丰度\n{receive_date}")

    # 整理数据格式
    rows = []
    for var in snv:
        gene = var.get("gene_symbol", "-")
        variant = (f'{var.get("hgvs_c", "-")}\n'
                   f'{var.get("hgvs_p", "-")}')
        if tissue_sample:
            tissue_sample_id = tissue_sample["sample_id"]
            tissue_freq = float(var[tissue_sample_id].replace("%", ""))/100 if re.search("%", str(var[tissue_sample_id])) else float(var[tissue_sample_id]) if isinstance(var[tissue_sample_id], str) and var[tissue_sample_id] else var[tissue_sample_id]
            tissue_freq = decimal_percen(tissue_freq) if tissue_freq else ""
        row = [gene, variant, tissue_freq]

        for sample in blood_sample:
            sample_id = sample["sample_id"]
            blood_freq = float(var[sample_id].replace("%", ""))/100 if re.search("%", str(var[sample_id])) else float(var[sample_id]) if isinstance(var[sample_id], str) and var[sample_id] else var[sample_id]
            blood_freq = decimal_percen(blood_freq) if blood_freq else ""
            row.append(blood_freq)
        rows.append(row)

    for var in sv:
        gene = var.get("gene_symbol", "-")
        variant = (f'{var.get("five_prime_gene")}-{var.get("three_prime_gene")}融合')
        if tissue_sample:
            tissue_sample_id = tissue_sample["sample_id"]
            tissue_freq = float(var[tissue_sample_id].replace("%", ""))/100 if re.search("%", str(var[tissue_sample_id])) else float(var[tissue_sample_id]) if isinstance(var[tissue_sample_id], str) and var[tissue_sample_id] else var[tissue_sample_id]
            tissue_freq = decimal_percen(tissue_freq) if tissue_freq else ""
        row = [gene, variant, tissue_freq]

        for sample in blood_sample:
            sample_id = sample["sample_id"]
            blood_freq = float(var[sample_id].replace("%", ""))/100 if re.search("%", str(var[sample_id])) else float(var[sample_id]) if isinstance(var[sample_id], str) and var[sample_id] else var[sample_id]
            blood_freq = decimal_percen(blood_freq) if blood_freq else ""
            row.append(blood_freq)
        rows.append(row)

    return {"headers": headers, "rows": rows}

def add_cell_to_row(row):
    """
    给 Word 表格行增加一个新的单元格
    """
    tc = OxmlElement("w:tc")
    row._tr.append(tc)
    return row.cells[-1]

def set_cell_text(cell, text, bold=False):
    """
    设置 Word 单元格内容
    """
    cell.text = ""
    paragraph = cell.paragraphs[0]
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run(str(text))
    run.bold = bold
    run.font.size = Pt(9)
    cell.vertical_alignment = (
        WD_CELL_VERTICAL_ALIGNMENT.CENTER
    )

def fill_table(table, table_data):
    """
    根据 headers + rows 填充 Word 表格
    """
    headers = table_data["headers"]
    rows = table_data["rows"]
    # ==================================================
    # 1. 删除模板中的数据行
    # ==================================================
    while len(table.rows) > 1:
        tr = table.rows[-1]._tr
        table._tbl.remove(tr)

    # ==================================================
    # 2. 处理表头
    # ==================================================
    header_row = table.rows[0]
    # 模板默认有3列：
    # 基因
    # 检测变异
    # 组织丰度
    # 根据 headers 动态增加后面的血检列
    # 先设置已有的3列
    for i in range(3):
        set_cell_text(
            header_row.cells[i],
            headers[i],
            bold=True
        )

    # 增加动态列
    for header in headers[3:]:
        cell = add_cell_to_row(header_row)
        set_cell_text(cell, header, bold=True)

    # ==================================================
    # 3. 添加数据
    # ==================================================
    for row_data in rows:
        row = table.add_row()
        # 前3列
        for i in range(3):
            set_cell_text(
                row.cells[i],
                row_data[i]
            )
        # 动态血检列
        for value in row_data[3:]:
            cell = add_cell_to_row(row)
            set_cell_text(
                cell,
                value
            )

def find_table(doc, marker):
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                if marker in cell.text:
                    return table
    raise ValueError(f"未找到包含标记 '{marker}' 的表格")