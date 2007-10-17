var get_jobs = new Ajax.PeriodicalUpdater('jobs_table', '/get_jobs',
  {
    method: 'get',
    frequency: 1,
    decay: 2
  });