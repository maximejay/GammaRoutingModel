source ~/python_venv/smash_new/bin/activate
nohup python3 exp_cance_ardeche.py --mode "init" --catchment "Cance" --scenario_key "s1" > log_init_cance_s1.txt &
nohup python3 exp_cance_ardeche.py --mode "smash" --catchment "Cance" --scenario_key "s1" > log_smash_cance_s1.txt &
nohup python3 exp_cance_ardeche.py --mode "smash_gamma" --catchment "Cance" --scenario_key "s1" > log_smash_gamma_cance_s1.txt &
