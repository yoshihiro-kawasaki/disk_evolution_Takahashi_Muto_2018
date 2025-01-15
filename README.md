# 原始惑星系円盤内のガスとダストの動径方向１次元計算 (Takahashi and Muto 2018)
Takahashi and Muto (2018) https://ui.adsabs.harvard.edu/abs/2018ApJ...865..102T/abstract を参照。
コードは、C++で書かれている。plotにpython(jupiter notebook)による結果の図示コードがある。

## 使い方
### コンパイル
Makefile内において、
```
PROBDIR = TM2018
```
の場合、ダストはストークス数をパラメータとする計算(Takahashi and Muto (2018))を行う。このとき、inputファイルは、input_TM2018を用いる。
```
PROBDIR = TM2018_2
```
の場合は、ダストサイズをパラメータとして計算を行う。Takahashi and Muto (2018) のFig. 16と17の計算を行う。
inputファイルは、input_TM2018_2を用いる。
```
PROBDIR = TM2018_dustgrowth
```
は、Takahashi and Muto (2018)の計算にsingle size approachによるダストサイズ進化(Sato et al. (2016), https://ui.adsabs.harvard.edu/abs/2016A%26A...589A..15S/abstract) を組み込んだものである。inputファイルは、input_TM2018_dustgrowthを用いる。
PROBDIRの選択後、makeでコンパイル。
```
$ make clean
$ make
```
`run`が生成される。

output用のディレクトリを作る。
```
$ mkdir output_dir_name
```

inputファイルの設定を行う。例えば、input/input_TM2018は以下である。
```
outdir = output/fiducial

# disk parameters
nr      = 150                    # number of radial cells
rmin    = 1.0e0                  # minimum (inner) radius [au]
rmax    = 1.0e4                  # maximum (outer) radius [au]
spacing = log                    # log or root

alpha_turb = 3.0e-4              # turbulent strength
c_wind     = 1.0e-4              # wind strength
St         = 1.0e-1              # Stokes number of dust
fdg        = 1.0e-2              # dust to gas mass ratio

# cloud parameters
nshells  = 50000                 # number of cloud shells
n_center = 2.555e5               # central number density of cloud [cm^-3]
Tc = 10.0                        # cloud temperature [K]
omega_c = 9.72234e-15            # cloud angular frequency [s^-1]
fc = 1.4                         # density factor to promote gravitational collapse takahashi et al 2013
Mc = 1.5                         # total cloud mass [solar mass]
Ms_init = 0.01                   # initial star mass [solar mass]

# time parameters
tend = 1.28e6                    # simulation end time [yr]
cfl  = 0.6                       # cfl number
delta_tout = 1.0e4               # delta output time [yr]
```
Takahashi and Muto (2018)では、特に乱流の強度`alpha_turb`、円盤風の強度`c_wind`、ダストのストークス数`St`をパラメータとして計算を行なっている。

計算の実行は、input_TM2018の場合
```
$ ./run input/input_TM2018
```
で行う。計算が開始すると、先ほど作成したoutput用のディレクトリ内に計算結果が出力される。

### 図示
plot内にpython(jupiter notebook)で書かれた簡単な図示するコード`plot.ipynb`がある。
inputファイルでoutputに設定したディレクトリのパスを設定する。
```
dir_name = "../output/fiducial/"
```
また、プロットしたい時間を４つ選択する
```
plot_time = [2.7e5, 6.0e5, 7.7e5, 1.27e6] # yr
```
後はコードを実行していくと面密度などが表示される。
他の物理量をプロットしたいときは、面密度のプロットを参考にする。