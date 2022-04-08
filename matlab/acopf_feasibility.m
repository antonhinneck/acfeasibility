solutions = csvread('pglib_opf_case118_ieee2022-02-19T14:01:41.506.txt');

allResults = zeros(2,0);
nsols = length(solutions(:, 1));

for k = 1:nsols

  sol_id = k;

  mpc = loadcase('case118');
  open_lines = zeros(0);
  ctr = 1;

  for i = 1:length(solutions(sol_id,:))
     mpc.branch(i,11) = solutions(sol_id, i);
     if solutions(sol_id, i) == 0
      open_lines(ctr) = i;
      ctr += 1;
     end
  end

  disp(ctr);

  [groups, isolated] = find_islands(mpc);
  find_islands(mpc);
  lab = build_lab(mpc);

  global points;
  points = zeros(length(isolated), 0);
  ## Function requires a global points variable.
  recursiveEnumPoints([], isolated, lab, points);
  disp(points);
  
  ctr = 1;
  successful_restorations = zeros(0,0);
  
  for i = 1:length(points(1,:))
    mpc1 = mpc;
    mpc1.branch(points(1, i), 11) = 1;
    mpc1.branch(points(2, i), 11) = 1;    
    [results, success] = runopf(mpc1);
    if success
      allResults(1,k) = 1;
      allResults(2,k) = results.f;
      successful_repair(:,ctr) = [i,j];
      ctr += 1;
    else
      allResults(1,k) = 0;
      allResults(2,k) = -1;
    end
  end
end

dlmwrite ("118_ac_feasibility.csv", allResults, ",")
