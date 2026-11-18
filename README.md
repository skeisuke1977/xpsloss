# xpsloss-starter

公開理論に基づく、XPSの誘電損失計算の独立実装に向けた着手版です。
作成日: 2026-11-19。

**現時点で動くのは、振動子で表した損失関数と、均質・等方な無限媒質のバルク非弾性散乱です。QUEELS-XPSの完全再現版ではありません。**

このパッケージはQUEELS/QUASESの公式製品ではなく、その実行ファイル・DLL・ソースコードを利用していません。公開されている数式から新規に実装しています。

## 現在の範囲

| 機能 | 状態 |
| --- | --- |
| JSON入力の検証 | 実装済み |
| Drude-Lindhard型振動子からELFを計算 | 実装済み |
| 運動量分散を含むバルクDIIMFP | 実装済み |
| CSV、JSONによる出力、計算条件の記録 | 実装済み |
| 解析解・和則・収束との照合 | テストあり |
| 複素誘電関数の因果的な再構成 | 未実装 |
| 表面、真空中の損失、静止正孔、干渉 | 未実装 |
| 深さ平均と自己無撞着な平均自由行程 | 未実装 |
| 多重散乱、スペクトル畳み込み、装置関数 | 未実装 |
| Quanty接続、CoO物性値、フィッティング | 未実装 |
| 元のQUEELS-XPSとの数値一致 | 未検証 |

ELFは Im[-1/epsilon]、DIIMFPは単位経路長・単位損失エネルギー当たりの非弾性散乱確率密度です。ELF単体を、測定されるXPSの損失確率と同一視しないでください。

## Windows PowerShellでの開始

Python 3.10以上を想定します。展開した `xpsloss-starter` フォルダーを開きます。

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m xpsloss capabilities --json
.\.venv\Scripts\python.exe -m xpsloss validate examples/synthetic_bulk.json
.\.venv\Scripts\python.exe -m xpsloss elf examples/synthetic_bulk.json --out runs/elf-demo
.\.venv\Scripts\python.exe -m xpsloss bulk examples/synthetic_bulk.json --out runs/bulk-demo
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

仮想環境のactivateは不要です。通常のインストールには依存パッケージ取得のため通信が必要ですが、計算処理そのものは外部APIを使いません。依存関係を事前に導入すればオフラインで計算できます。

## macOS / Linux

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e .
.venv/bin/python -m xpsloss capabilities --json
.venv/bin/python -m xpsloss bulk examples/synthetic_bulk.json --out runs/bulk-demo
.venv/bin/python -m unittest discover -s tests -v
```

**例題は合成材料です。Al、CoO、Siなどの実測・文献データとして利用しないでください。** 既存の出力ディレクトリには上書きしません。再計算するときは新しい出力名を指定します。

## AIエージェントとの共同開発

最初に `AGENTS.md`、`docs/PROJECT_PLAN.md`、`docs/EQUATIONS.md` を渡してください。`CLAUDE.md`にも同じ資料への案内があります。

次に依頼する仕事は `docs/NEXT_TASK.md` に書いてあります。バルクの実装を足場として、まず1997年原著の数式と計算手順の対応を確定します。ソフトの機能が増えたことと、物理モデルの妥当性が検証されたことは、別々に記録します。

## 検証と資料

- `docs/VALIDATION.md`: この配布物で実行した検証と未検証事項。
- `docs/EQUATIONS.md`: 数式、単位、検証の根拠。
- `docs/PROJECT_PLAN.md`: 完全モデルとQuanty接続へ進むための設計。
- `docs/REFERENCES.md`: 一次資料と用途。
- `NOTICES.md`: 独立実装、参照資料、公開時の整理。

WindowsとmacOSでの実行結果は、この初版ではまだ確認していません。検証済みの実行環境は `docs/VALIDATION.md` を参照してください。
