# AutoReport_v5 新人上手指南

## 1. 代码库整体结构

- `main.py`：外部调用入口，负责参数解析并调用核心处理函数。  
- `AutoReport_Base_Code_v1_0_0.py`：主流程编排器，完成 JSON 读取、模板匹配、分模块数据处理、Word 报告渲染。  
- `bin/`：报告内容构建模块，按业务域拆分（样本信息、QC、变异、MSI/TMB/HRD/TME、药物与临床试验、参考文献等）。  
- `libs/`：底层工具函数和规则处理（变异格式化、证据规则、配置读取、排序规则、检测结果汇总等）。  
- `merge_product/`：联合产品/合并订单的特化逻辑（例如 BRCA+HRD、DNA+RNAseq 等）。  
- `image/`：报告中使用的固定图像素材。  
- `Auto_test/JSON_SET/`：自动化测试样例数据（按产品线组织）。

## 2. 关键执行链路（建议先理解）

1. `main.py` 调用 `AutoReport_Base_Code_v1_0_0.get_data(...)`。  
2. `get_data` 读取 `<outfile>/<json_name>.json` 作为输入。  
3. 如果命中 `merge_order` 或特定 `relation_order_code` 场景，优先走 `merge_product/` 的合并逻辑。  
4. 常规流程下依次执行：
   - 产品别名与模块类型归一化；
   - 证据裁剪（MSS/TP53/KRAS 等按规则过滤）；
   - 模板选择（`getTemplate.choose_template`）；
   - 各业务模块取数并写入 `data`；
   - 固定/动态参考文献组装；
   - 最终渲染并输出报告文件（Word/JSON 视配置）。

## 3. 新人必须优先掌握的 8 个点

1. **输入 JSON 契约**：`sample_info`、`snvindel/cnv/sv`、`merge_order` 等字段结构。  
2. **模板匹配规则**：`bin/getTemplate.py` 与 `config` 下映射关系如何决定 `report_name`。  
3. **变异处理主干**：`bin/getVar.py` 与 `libs/processThe*` 系列的分型/注释/证据级别处理。  
4. **证据删减策略**：`bin/remove*.py` 中按产品、癌种、公司定制的业务规则。  
5. **联合产品逻辑**：`merge_product/` 中“多样本合单报告”触发条件与字段拼接方式。  
6. **参考文献来源**：`bin/getRefence.py` 的 fixed/dynamic 生成机制，便于审阅报告可追溯性。  
7. **测试样例组织**：`Auto_test/JSON_SET/` 对应的产品线和场景覆盖。  
8. **编码与字符集历史包袱**：文件存在 `gbk` 声明与中文注释，改动时要避免编码损坏。

## 4. 推荐学习路径（2 周）

### 第 1-2 天：跑通最小链路
- 阅读 `main.py`、`AutoReport_Base_Code_v1_0_0.py` 的入口与主函数。  
- 随机选一个 `Auto_test/JSON_SET/` 场景，手动追踪数据在 `data` 中的落点。

### 第 3-5 天：吃透模板与变异
- 深读 `bin/getTemplate.py`、`bin/getVar.py`。  
- 对照 `libs/processTheSnvindel.py`、`libs/processTheCNV.py`、`libs/processTheSV.py` 看不同变异类型的处理差异。

### 第 6-8 天：理解临床解释与证据裁剪
- 阅读 `bin/getApprovalDrug.py`、`bin/getClinicTrial.py`、`bin/removeEvi.py` 及相关模块。  
- 汇总“哪些条件会删证据”的规则清单，形成个人笔记。

### 第 9-10 天：学习合并产品
- 阅读 `merge_product/` 目录下已有案例。  
- 重点关注触发条件：产品组合、医院（company/origin_company）、关联订单字段。

### 第 11-14 天：小步改造 + 回归
- 先做“低风险改动”（日志、注释、配置键兼容）。  
- 每次改动后用至少 2-3 组 `Auto_test/JSON_SET/` 样例回归，确认报告结构与关键信息不回退。

## 5. 实操建议（避免踩坑）

- 先**加样例再改逻辑**：新增规则前先准备最小 JSON 复现。  
- 业务规则尽量模块内聚：产品/医院定制逻辑优先放对应 `remove*` 或 `merge_product`，不要散落主流程。  
- 改字段名前先做兼容：此仓库与上游系统有耦合，字段变更要保留兜底。  
- 关注“报告正确性优先于代码优雅”：这是生产型报告引擎，宁可保守改造也别大面积重构。

## 6. 给新人的第一批阅读清单（按优先级）

1. `main.py`  
2. `AutoReport_Base_Code_v1_0_0.py`  
3. `bin/getTemplate.py`  
4. `bin/getVar.py`  
5. `bin/getRefence.py`  
6. `libs/getConfig.py`  
7. `merge_product/*.py`（按当前负责产品线选择）
