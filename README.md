# nmr-1h-duplicate-check

已经检查并补齐了功能。现在这个 skill 支持：

1. 下载后让 Codex 识别并安装。
2. 读取 `.mnova` 中批量保存的 ¹H NMR 引用数据。
3. 对同一个 `.mnova` 内部的谱图做两两查重。
4. 使用固定的 `--history-dir` 后，新上传文件会自动和之前所有已登记文件一起比较。
5. 输出 `NEW_NEW`、`NEW_HISTORY`、`HISTORY_HISTORY` 三种比较关系，方便确认新文件确实和历史文件做了比较。
6. 当前 skill 是读取 `.mnova` 中记录的 **核磁共振氢谱原始数据路径**，不是直接解析 `.mnova` 内嵌的谱图。因此 `.mnova` 换到另一台电脑时，原始数据目录也必须能访问；否则这条路径无法读取。
