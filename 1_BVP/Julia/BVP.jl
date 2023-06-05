function thomas(a,b,c,d)
    last_x_index = length(a)  # 要素数を算出
    # Thomas法の算出する係数
    g = zeros(Float64, last_x_index)
    s = zeros(Float64, last_x_index)
    u = zeros(Float64, last_x_index)
    # 係数を算出していく
    for i = first_x_index:last_x_index
        if i==1
            g[i] = b[i]
            s[i] = d[i]
        else
            g[i] = b[i] - a[i]*c[i-1] / g[i-1]
            s[i] = d[i] - a[i]*s[i-1] / g[i-1]
        end
    end
    # 一番最後の方程式の解を求める
    u[last_x_index] = s[last_x_index] / g[last_x_index]
    # 解を求める 
    for j = reverse(first_x_index:last_x_index)  # 逆からループする
        if j==last_x_index
            continue
        end
        u[j] = (s[j] - c[j]*u[j+1]) / g[j] 
    end
    return u
end

n = 20 # mesh number for x
Δx = 1/n
# 元のコードは、数式[0,n] -> [1,n+1]の格子になっていた。
# おそらく、Fortranが1インデックスから始まるためにそれに合わせたのだと思われる。
# しかし、ややこしいので、[0,n] の格子を使う。境界条件の 0,nは常に0なので省いて、[1,n-1]の配列を常に考える。
first_x_index = 1  # il
last_x_index = n-1  # i

# n個の要素で0で初期化したベクトルを作成。Thomas法の係数部分に当たる。別にこういうふうに書く必要はあまりない。
a = zeros(Float64, last_x_index)
b = zeros(Float64, last_x_index)
c = zeros(Float64, last_x_index)
d = zeros(Float64, last_x_index)
# 係数行列と照らし合わせて、Thomas法の係数に当たるものを入れていく。
# i は行の番号。
for i = first_x_index:last_x_index  # i = 1,2,...,n-1
    a[i] = 1
    b[i] = Δx^2 - 2
    c[i] = 1
    d[i] = -i * Δx^3
end
# トーマス法で解く
u = thomas(a,b,c,d)

# 解析解を作る
u_analytical = Vector{Float64}()
x_array = Vector{Float64}() 
for i = first_x_index:last_x_index
    x = Δx * i
    push!(x_array, x)
    u_ = sin(x)/sin(1) - x
    push!(u_analytical, u_)  # 空の配列 u_analytical に u_ を追加する。
end

# A / B は Aに右からBの逆行列を掛ける操作になる。
# A ./ B は要素同士の割り算。　ドット演算はドットをつける。
# 元のコードと同じエラーを算出
err = (u - u_analytical) ./ u_analytical * 100

#=
# ファイルを書き出す。 適当なら以下でOK
open("output.dat", "w") do out
    Base.print_array(out, hcat(x_array[:], u[:], u_analytical[:], err[:]))
end
=#

# ちゃんと書き出し方を指定したい場合は以下
using Printf
open( "output.dat", "w" ) do out
    [@printf(out, "%7.7e %7.7e %7.7e %7.7e\n", x_array[k], u[k], u_analytical[k], err[k]) for k=1:length(x_array[:])]
end


#=ライブラリがない場合入れる。
# ターミナルで、Macなら、
brew install ffmpeg

juliaのREPLから、
using Pkg
Pkg.add("FFMPEG")
Pkg.add("ImageMagick")
Pkg.build("ImageMagick")
Pkg.add("Plots")
Pkg.update()
=#
# プロットする。
using Plots
p = plot()
scatter!(x_array, u, alpha=0.5, line_width=5, label="u")
plot!(x_array, u_analytical, linestyle=:dash, label="u_analytical")
scatter!(x_array, err, alpha=0.5, line_width=5, label="err")
savefig(p, "BVP.png")