function pool = start_optional_thread_pool(cfg)
pool = gcp('nocreate');
if cfg.use_thread_pool
    if isempty(pool)
        pool = parpool('Threads', cfg.thread_pool_workers);
    elseif ~isa(pool, 'parallel.ThreadPool')
        warning('An existing non-thread pool is active; leaving it unchanged.');
    end
end
end
