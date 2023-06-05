using Plots

function heat_e(boundary_condition = "Dirichlet-0")
    # parameter
    n_mesh = 20  # メッシュ数
    n_step = 250  # ステップ数
    α = 1  # 熱伝達率
    Δt = 0.001  # 時間ステップ
    Δx = 1/(n_mesh-1)  # 格子幅
    r = α * Δt / Δx^2  # 元の数式を参照
    mesh_half_point = Int(trunc((n_mesh + 1)/ 2))  # 初期条件等を考えるとき、奇数の場合は最初の方を多めにする。

    # initial condition
    u = zeros(Float64, n_mesh)
    x_array = []
    for j = 1:n_mesh
        x = (j - 1) / (n_mesh - 1)
        push!(x_array, x)
        if j <= mesh_half_point
            u[j] = x
        else
            u[j] = 1 - x
        end
    end


    # main loop
    u_save = []
    # initialization
    uⁿ⁺¹ = zeros(Float64, n_mesh)

    for n = 1:n_step

        # boundary condition
        u[1] = 0.0
        if boundary_condition == "Dirichlet-0"
            u[n_mesh] = 0.0  # "Dirichlet-0"
        elseif boundary_condition == "Neumann-0"
            u[n_mesh] = 2r * u[n_mesh-1] + (1-2r)*u[n_mesh]  # Neumann-0
        elseif boundary_condition == "Neumann-1"
            u[n_mesh] = 2r * u[n_mesh-1] + (1-2r)*u[n_mesh] + 2r*Δx  # Neumann-1
        else
            throw(DomainError(boun, "boundary_condition is not correct"))
        end
        
        # save
        push!(u_save, u)

        for j = 2:n_mesh-1
            uⁿ⁺¹[j] = r*u[j-1] + (1-2r)u[j] + r*u[j+1]
        end

        # next step
        u = copy(uⁿ⁺¹)
        if n==n_step
            push!(u_save, u)
        end

        if abs(u[mesh_half_point]) >= 10000
            print("now is $(n) steps of $(n_step) steps, the calculation is diverge!")
        end
    end

    mkpath("anime")
    anim = Animation()
    for n in 1:n_step
        plt = plot(xlims=(0,1), ylims=(0,1))
        plot!(x_array, u_save[n], label="u")
        frame(anim, plt)
    end
    gif(anim, "anime/heat_$(boundary_condition).gif", fps=20)
end

# 第一引数にboundary_conditionが入る。ターミナルで以下を実行。
# julia heat_e.jl Dirichlet-0
# julia heat_e.jl Neumann-0
# julia heat_e.jl Neumann-1
if abspath(PROGRAM_FILE) == @__FILE__  # if __name__=="__main__"
    boundary_condition = ARGS[1]
    heat_e(boundary_condition)
end