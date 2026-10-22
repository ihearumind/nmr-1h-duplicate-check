# nmr-1h-duplicate-check

下载后使用 Codex 识别并安装本 skill。输入为 MestReNova（`.mnova`）中批量保存的 ¹H NMR 谱；Codex 会读取其引用的 Bruker 数据并直接做文件内两两查重。使用固定的 `--history-dir` 运行时，新上传文件会自动与之前所有已登记文件一起比较。

## 安装

把 `nmr-1h-duplicate-check` 文件夹复制到：

```text
$CODEX_HOME/skills/
```

依赖：

```bash
python -m pip install numpy scipy openpyxl
```

## 使用

先从 `.mnova` 提取 Bruker 数据路径（原始数据目录需要仍可访问；提取后请核对 inventory 数量）：

```bash
python scripts/extract_mnova_paths.py \
  --mnova "A.mnova" "B.mnova" \
  --search-root "F:/path/to/NMR" \
  --out inventory.json
```

再运行区间查重：

```bash
python scripts/nmr_range_duplicate.py \
  --inventory inventory.json \
  --history-dir ./nmr_history \
  --out nmr_results \
  --xlsx
```

以后继续使用同一个 `--history-dir`，新文件会与历史库中的全部谱图一起重新比较。

## 输出

高置信重复、疑似重复、全部两两比较、区间证据、数据集 QC、排除项，以及按“新文件/历史文件”拆分的比较范围统计。

高相似 ¹H 谱不等于同一化合物；结构判断仍需结合 ¹³C、2D NMR、MS 和样品记录。
