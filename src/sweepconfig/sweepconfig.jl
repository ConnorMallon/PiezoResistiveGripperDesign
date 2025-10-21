
function save_job_dicts(sweep_params,fixed_params,job_id)
    # save array of dictionaries for sweeping
    dicts = dict_list(sweep_params); println("Total dictionaries made: ", length(dicts))
    dict_vector = tmpsave(dicts,projectdir("_tmp"))
    dict_vector_dict = Dict("dict_vector" => dict_vector)
    x = "dict_vector_$(job_id)"
    dict_vector_sname = savename((@dict x),"jld2")
    DrWatson.save(datadir("results", "dict_vector", dict_vector_sname), dict_vector_dict)

    # save array of fixed params
    y = "fixed_dict_vector_$(job_id)"*".jld2"
    fixed_params
    DrWatson.save(datadir("results", "dict_vector", y), fixed_params)
end

function pbs_sweeper(sweep_params,fixed_params)
    job_id = rand(1:100000000)
    save_job_dicts(sweep_params,fixed_params,job_id)

    # setup job submission parametes
    job_name = sweep_params["model_name"][1]
    array_length =  prod(length(val) for val in values(sweep_params))
    julia_bin = Base.julia_cmd()[1]
    #julia_bin = "/home/mallon2/.julia/juliaup/julia-1.11.3+0.x64.linux.gnu/bin/julia"
    #julia_bin = "/home/565/cm8825/.julia/juliaup/julia-1.11.2+0.x64.linux.gnu/bin/julia"

    std_root_path1 = datadir("output_files")
    std_root_path2 = datadir("output_files")
    runcase_path = srcdir("sweepconfig/runcase.jl")
    runcase_path2 = srcdir("sweepconfig/runcase2.jl")

    project_dir = projectdir()

    job_array_limit = 10

for start_idx in 1:job_array_limit:array_length

    end_idx = min(start_idx + job_array_limit - 1, array_length)

    @show job_name

    runcase2_path = srcdir("sweepconfig/runcase2.jl")

    ncpus = 1
    mem = 4*4 # GB
    time = 24 # hours

    mem_string = "$(mem)G"
    hours_string = "$(time):00:00"

    job_cost = array_length*ncpus*(mem/4)*time*2
    println("cost of job (SU): ",job_cost) # 2 for normal queue 

    if job_cost > 40000
        @show job_cost
        error("job cost too high, exiting")
    end




    bash_script = 
    """
    #!/bin/bash
    #
    #PBS -N $job_name
    #
    #PBS -o /dev/null
    #PBS -e /dev/null
    #
    #PBS -l walltime=24:00:00
    #
    #PBS -l mem=$mem_string
    #PBS -l ncpus=$ncpus
    #
    #PBS -J $start_idx-$end_idx
    #
    #PBS -r y

    output_file="$(std_root_path2)/file-\${PBS_JOBID}.ext1"
    error_file="$(std_root_path2)/file-\${PBS_JOBID}.ext2"
    exec > "\$output_file" 2> "\$error_file"

    cd $project_dir
    pwd
    $julia_bin $runcase_path $job_id \${PBS_ARRAY_INDEX}

    """

    bash_script_file = tempname()
    write(bash_script_file, bash_script)
    read(run(`qsub $(string(bash_script_file))`), String)
    rm(bash_script_file)
end

end
