# 独立実装と資料

この着手版は、公開された数式を参照して新規に作ったコードと開発資料から成る。QUEELS-XPSの公式版・承認版・互換保証版ではない。

元製品の実行ファイル、DLL、ソース、図版、論文PDFは含めていない。論文は共有の論文庫のキー（例 `Tougaard2022`、`Tougaard2011`）で参照し、PDFやページ画像はリポジトリに入れない（`.gitignore` で除外）。

## 同梱データの区分

| 場所 | 内容 | 由来と条件 |
| --- | --- | --- |
| `examples/` | 合成パラメータ | 実材料に由来しない。Al、CoO、Si等の実測・文献値として使わない |
| `benchmarks/sources/zenodo6024064_v1.0/` | ELFデータ集の原本3ファイル（無改変、バイト列とSHA-256を保持） | N. Pauly, F. Yubero, S. Tougaard, “ELF dielectric functions for various materials”, Version 1.0, Zenodo, 2022-02-09, https://doi.org/10.5281/zenodo.6024064 。[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| `benchmarks/inputs/*elfdb*` | 上記データ集からSi・Au・Alの材料パラメータを抜粋し、来歴を付けてJSONに転記したもの。強度の再規格化はしていない | 同上（CC BY 4.0）。各記録に原本の行番号とSHA-256を記録 |
| `benchmarks/inputs/si_queels_manual_v2.1.json` | QUEELS-XPS v2.1 manual p. 26 表1とp. 27の設定値を、ページ番号付きの事実として記録 | manualは S. Tougaard, F. Yubero の著作物（copyright 2001–2011、権利留保）。本文・図は複製していない |

ELFデータ集の値は文献から収集された物性入力であり、本パッケージがその値の妥当性を検証したものではない。実材料の入力を同梱していても、QUEELS-XPSとの数値一致は未検証である（[VALIDATION.md](docs/VALIDATION.md)）。

元製品のページとmanualには非商用利用無料・権利留保の記載がある。ダウンロード無料という事実を、コード転用や再配布の許諾と同一視しない。外部データの追加時は、そのデータに適用される条件を個別に記録する。

## 公開前に決めること

コードの配布ライセンス（`LICENSE`）は未確定。確定時には、`benchmarks/sources/**` の原本と `benchmarks/inputs/*elfdb*` の派生データをコードライセンスの適用範囲から除外し、CC BY 4.0 が適用される旨と上記の帰属表示を明記する。

## 参照ソフトの出力（決定済み）

QUEELS-XPS で計算した参照出力（`work/reference_runs/`）は**公開しない**（利用者の決定、2026-11-19）。gitの管理外に置き、リポジトリにも配布物にも含めない。文書では、QUEELS-XPS とその理論の原著（Simonsen, Yubero, Tougaard 1997 など）を出典として引用し、比較の条件と結果（許容幅に対する差の指標など）を書いてよい。出力ファイルそのものや、数値表をまとめて転載することはしない。
