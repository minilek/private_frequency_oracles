class Params:
    

    def __init__(self, input_file=None, output_dir=None, num_users=1000, universe_size=8192,
                 num_trials=400, epsilon=5.0,  debug=False, include_ss=True,
                 include_rappor=True, include_pirappor=True, include_rr=True,
                 include_hr=True, include_vhr=False, include_rhr=True, include_pg=True,
                 include_hpg=True, run_variance_experiments=True,
                 run_timing_experiments=True, run_mse_experiments=True,
                 variance_repetitions=10, var_epsilon_min=1.1, var_epsilon_max=5.5,
                 num_var_epsilons=45, universe_min=16, universe_max=8192, hpgq=5,
                 num_universes=100):
        self.num_users = num_users
        self.universe_size = universe_size
        self.num_trials = num_trials
        self.input_file = input_file
        self.output_dir = output_dir
        self.debug = debug
        self.include_ss = include_ss
        self.include_rappor = include_rappor
        self.include_pirappor = include_pirappor
        self.include_rr = include_rr
        self.include_hr = include_hr
        self.include_vhr = include_vhr
        self.include_rhr = include_rhr
        self.include_pg = include_pg
        self.include_hpg = include_hpg
        self.run_variance_experiments = run_variance_experiments
        self.run_timing_experiments = run_timing_experiments
        self.run_mse_experiments = run_mse_experiments
        if run_mse_experiments or run_timing_experiments:
            self.epsilon = epsilon
        if run_variance_experiments:
            self.variance_repetitions = variance_repetitions
            self.var_epsilon_min = var_epsilon_min
            self.var_epsilon_max = var_epsilon_max
            self.num_var_epsilons = num_var_epsilons
        if run_timing_experiments:
            self.universe_min = universe_min
            self.universe_max = universe_max
            self.num_universes = num_universes
        self.hpgq = hpgq

    def __str__(self):
        cmd = './run_experiment'
        cmd += ' --num-users=' + str(self.num_users)
        if self.input_file!=None: cmd += ' --input-file=' + self.input_file
        if self.output_dir!=None: cmd += ' --output-dir=' + self.output_dir
        if self.run_mse_experiments or self.run_timing_experiments:
            cmd += ' --epsilon=' + str(self.epsilon)
        if self.run_mse_experiments or self.run_variance_experiments:
            cmd += ' --universe-size=' + str(self.universe_size)
        if self.run_mse_experiments: cmd += ' --num-trials=' + str(self.num_trials)
        if self.debug: cmd += ' --debug'
        if self.include_ss: cmd += ' --include-ss'
        if self.include_rappor: cmd += ' --include-rappor'
        if self.include_pirappor: cmd += ' --include-pirappor'
        if self.include_rr: cmd += ' --include-rr'
        if self.include_hr: cmd += ' --include-hr'
        if self.include_vhr: cmd += ' --include-vhr'
        if self.include_rhr: cmd += ' --include-rhr'
        if self.include_pg: cmd += ' --include-pg'
        if self.include_hpg: cmd += ' --include-hpg --hpgq=' + str(self.hpgq)
        if self.run_variance_experiments:
            cmd += ' --run-variance-experiments'
            cmd += ' --variance-repetitions=' + str(self.variance_repetitions)
            cmd += ' --var-epsilon-min=' + str(self.var_epsilon_min)
            cmd += ' --var-epsilon-max=' + str(self.var_epsilon_max)
            cmd += ' --num-var-epsilons=' + str(self.num_var_epsilons)
        if self.run_timing_experiments:
            cmd += ' --run-timing-experiments'
            cmd += ' --universe-min=' + str(self.universe_min)
            cmd += ' --universe-max=' + str(self.universe_max)
            cmd += ' --num-universes=' + str(self.num_universes)
        if self.run_mse_experiments: cmd += ' --run-mse-experiments'
        return cmd
            
    @staticmethod
    def string_to_param(s):
        params = s.split(' ')[1:]
        p = Params(include_ss=False, include_rappor=False, include_pirappor=False,
                   include_rr=False, include_hr=False, include_rhr=False,
                   include_pg=False, include_hpg=False, run_variance_experiments=False,
                   run_timing_experiments=False, run_mse_experiments=False)
        for param in params:
            l = param.split('=')
            option = l[0][2:].replace('-', '_') # '2:' removes the -- at beginning of option name
            if len(l) == 1:
                p.__dict__[option] = True
            elif l[1].isdigit():
                p.__dict__[option] = int(l[1])
            else:
                try:
                    f = float(l[1])
                    p.__dict__[option] = f
                except ValueError:
                    p.__dict__[option] = l[1]
        return p

                    
