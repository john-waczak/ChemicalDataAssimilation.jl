using Pkg
Pkg.activate(".")
using ArgParse


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--mechanism_path"
            help = "Path to mechanism `.fac` file specifying the chemical mechanism to be used."
            arg_type = String
            default = "mechanism-files/extracted/alkanes/methane.fac"
        "--model_name", "-n"
            help = "Name for the resulting model used in output paths"
            arg_type = String
            default = "methane"
        "--time_step", "-t"
            help = "The time step used during integration of mechanism (in minutes)."
            arg_type = Float64
            default = 15.0
        "--use_updated_photolysis"
            help = "Whether or not to use updated photolysis"
            arg_type = Bool
            default = true
        "--restart"
            help = "Whether or not to restart 4d_var from previous fitresult"
            action = :store_true
        "--try_solve"
            help = "Whether or not to precompile solvers by calling once."
            action = :store_true
        "--use_background_cov"
            help = "Whether or not to use background covariance matrix in loss"
            action = :store_true
        "--fudge_fac", "-f"
            help = "A fudge factor for manipulating scale of measurement uncertainties"
            arg_type = Float64
            default = 0.5
        "--epsilon", "-e"
            help = "Estimated background uncertainty for diagonal of B matrix, i.e. uncertainty in initial condition"
            arg_type = Float64
            default = 0.5
    end

    parsed_args = parse_args(ARGS, s; as_symbols=true)
    return parsed_args
end



function create_slurm_scripts(parsed_args; n_tasks=8)
    preamble = """
    #!/bin/bash

    #SBATCH     --job-name=$(parsed_args[:model_name])
    #SBATCH     --output=$(parsed_args[:model_name]).out
    #SBATCH     --error=$(parsed_args[:model_name]).err
    #SBATCH     --nodes=1
    #SBATCH     --ntasks=1
    #SBATCH     --cpus-per-task=$(n_tasks)   # number of threads for multi-threading
    #SBATCH     --time=2-00:00:00  # 2 day max
    #SBATCH     --mem=12G
    #SBATCH     --mail-type=ALL
    #SBATCH     --mail-user=jxw190004@utdallas.edu
    #SBATCH     --partition=normal

    """


    step1 = """
        julia --threads \$SLURM_CPUS_PER_TASK --project=. 1__build_simulation.jl --mechanism_path $(parsed_args[:mechanism_path])  --model_name $(parsed_args[:model_name]) --time_step $(parsed_args[:time_step]) --use_updated_photolysis $(parsed_args[:use_updated_photolysis])
        """

    step2 = """
        julia --threads \$SLURM_CPUS_PER_TASK --project=. 2__run_4dvar.jl --mechanism_path $(parsed_args[:mechanism_path])  --model_name $(parsed_args[:model_name]) --time_step $(parsed_args[:time_step]) --restart $(parsed_args[:restart]) --try_solve $(parsed_args[:try_solve]) --use_background_cov $(parsed_args[:use_background_cov]) --fudge_fac $(parsed_args[:fudge_fac]) --epsilon $(parsed_args[:epsilon])
        """

    step2b = """
        julia --threads \$SLURM_CPUS_PER_TASK --project=. 2b__visualize_results.jl --model_name $(parsed_args[:model_name]) --time_step $(parsed_args[:time_step]) --restart $(parsed_args[:restart]) --fudge_fac $(parsed_args[:fudge_fac])
        """

    step3 = """
        julia --threads \$SLURM_CPUS_PER_TASK --project=. 3__run_ekf.jl --mechanism_path $(parsed_args[:mechanism_path])  --model_name $(parsed_args[:model_name]) --time_step $(parsed_args[:time_step]) --try_solve $(parsed_args[:try_solve]) --fudge_fac $(parsed_args[:fudge_fac]) --epsilon $(parsed_args[:epsilon])
        """

    step3b = """
        julia --threads \$SLURM_CPUS_PER_TASK --project=. 3b__visualize_results.jl --model_name $(parsed_args[:model_name])
        """


    open("1__$(parsed_args[:model_name]).slurm", "w") do f
        println(f, preamble)
        println(f, step1)
    end

    open("2__$(parsed_args[:model_name]).slurm", "w") do f
        println(f, preamble)
        println(f, step2)
    end

    open("2b__$(parsed_args[:model_name]).slurm", "w") do f
        println(f, preamble)
        println(f, step2b)
    end

    open("3__$(parsed_args[:model_name]).slurm", "w") do f
        println(f, preamble)
        println(f, step3)
    end

    open("3b__$(parsed_args[:model_name]).slurm", "w") do f
        println(f, preamble)
        println(f, step3b)
    end


    final_script = """
    #!/bin/bash

    RES=\$(sbatch 1__$(parsed_args[:model_name]).slurm)
    RES2=\$(sbatch --dependency=afterany:\${RES##* }  2__$(parsed_args[:model_name]).slurm)
    sbatch --dependency=afterany:\${RES2##* }  2b__$(parsed_args[:model_name]).slurm)
    RES3=\$(sbatch --dependency=afterany:\${RES2##* }  3__$(parsed_args[:model_name]).slurm)
    sbatch --dependency=afterany:\${RES3##* }  3b__$(parsed_args[:model_name]).slurm)

    """
    open("submit_jobs__$(parsed_args[:model_name]).sh", "w") do f
        println(f, final_script)
    end
end


parsed_args = parse_commandline()
create_slurm_scripts(parsed_args)


