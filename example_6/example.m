function out = model
%
% meshsweep_v4.m
%
% Model exported on May 19 2026, 22:06 by COMSOL 5.6.0.280.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');

model.study.create('std1');
model.study('std1').create('eig', 'Eigenfrequency');
model.study('std1').feature('eig').set('conrad', '1');
model.study('std1').feature('eig').set('solnum', 'auto');
model.study('std1').feature('eig').set('notsolnum', 'auto');
model.study('std1').feature('eig').activate('ewfd', true);

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('ewfd').prop('ShapeProperty').set('order_electricfield', 1);

model.param.set('ax', '3[mm]');
model.param.descr('ax', '');
model.param.set('ay', '3[mm]');
model.param.descr('ay', '');
model.param.set('az', '4.5[mm]');
model.param.descr('az', '');
model.param.set('d_pcb', '0.75[mm]');
model.param.descr('d_pcb', '');
model.param.set('d_metal', '1.6[mm]');
model.param.descr('d_metal', '');
model.param.set('thick_metal', '0.2[mm]');
model.param.descr('thick_metal', '');
model.param.set('R_metal', '0.4[mm]');
model.param.descr('R_metal', '');
model.param.set('az_total', 'az+2*d_pcb');
model.param.descr('az_total', '');
model.param.set('ax', '3[mm]');
model.param.descr('ax', '');
model.param.set('ay', '3[mm]');
model.param.descr('ay', '');
model.param.set('az', '4.5[mm]');
model.param.descr('az', '');
model.param.set('d_pcb', '0.75[mm]');
model.param.descr('d_pcb', '');
model.param.set('d_metal', '1.6[mm]');
model.param.descr('d_metal', '');
model.param.set('thick_metal', '0.2[mm]');
model.param.descr('thick_metal', '');
model.param.set('R_metal', '0.4[mm]');
model.param.descr('R_metal', '');
model.param.set('az_total', 'az+2*d_pcb');
model.param.descr('az_total', '');

model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').set('unite', true);
model.component('comp1').geom('geom1').feature('wp1').geom.create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('size', {'ax' 'ay'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.run('r1');
model.component('comp1').geom('geom1').feature('wp1').geom.create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('c1').set('r', 'R_metal');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('c1').set('pos', {'d_metal/2' 'd_metal/2'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('c1');
model.component('comp1').geom('geom1').feature('wp1').geom.run('c1');
model.component('comp1').geom('geom1').feature('wp1').geom.create('arr1', 'Array');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('arr1').set('fullsize', [2 2]);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('arr1').set('displ', {'-d_metal' '-d_metal'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('c1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('arr1').selection('input').set({'c1'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('arr1');
model.component('comp1').geom('geom1').feature('wp1').geom.run('arr1');
model.component('comp1').geom('geom1').feature('wp1').geom.create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('size', {'d_metal+thick_metal*2' '1'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').setIndex('size', 'd_metal+thick_metal*2', 1);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r2').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.run('r2');
model.component('comp1').geom('geom1').feature('wp1').geom.create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('size', {'d_metal-thick_metal*2' '1'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').setIndex('size', 'd_metal-thick_metal*2', 1);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('r3').set('base', 'center');
model.component('comp1').geom('geom1').feature('wp1').geom.run('r3');
model.component('comp1').geom('geom1').feature('wp1').geom.create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('dif1').selection('input').set({'r2'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('dif1').selection('input2').set({'r3'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('dif1');
model.component('comp1').geom('geom1').feature('wp1').geom.create('dif2', 'Difference');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('dif2').selection('input').set({'dif1'});
model.component('comp1').geom('geom1').feature('wp1').geom.feature('dif2').selection('input2').set({'arr1'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('dif2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('dif2').set('keep', true);
model.component('comp1').geom('geom1').feature('wp1').geom.run('dif2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('dif2').set('keep', false);
model.component('comp1').geom('geom1').feature('wp1').geom.run('dif2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.duplicate('c2', 'c1');
model.component('comp1').geom('geom1').feature('wp1').geom.feature.duplicate('arr2', 'arr1');
model.component('comp1').geom('geom1').feature('wp1').geom.run('c2');
model.component('comp1').geom('geom1').feature('wp1').geom.run('c2');
model.component('comp1').geom('geom1').feature('wp1').geom.run('c2');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('arr2').selection('input').set({'c2'});
model.component('comp1').geom('geom1').feature('wp1').geom.run('arr2');
model.component('comp1').geom('geom1').run('wp1');
model.component('comp1').geom('geom1').feature.create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 0, 0);
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 'd_pcb', 1);
model.component('comp1').geom('geom1').feature('ext1').set('displ', {'0' '0'; '0' '0'});
model.component('comp1').geom('geom1').feature('ext1').set('scale', {'1' '1'; '1' '1'});
model.component('comp1').geom('geom1').feature('ext1').set('twist', {'0' '0'});
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 'd_pcb+az', 2);
model.component('comp1').geom('geom1').feature('ext1').set('displ', {'0' '0'; '0' '0'; '0' '0'});
model.component('comp1').geom('geom1').feature('ext1').set('scale', {'1' '1'; '1' '1'; '1' '1'});
model.component('comp1').geom('geom1').feature('ext1').set('twist', {'0' '0' '0'});
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 'd_pcb+az+d_pcb', 3);
model.component('comp1').geom('geom1').feature('ext1').set('displ', {'0' '0'; '0' '0'; '0' '0'; '0' '0'});
model.component('comp1').geom('geom1').feature('ext1').set('scale', {'1' '1'; '1' '1'; '1' '1'; '1' '1'});
model.component('comp1').geom('geom1').feature('ext1').set('twist', {'0' '0' '0' '0'});
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 0, 0);
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 'd_pcb', 1);
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 'd_pcb+az', 2);
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 'd_pcb+az+d_pcb', 3);
model.component('comp1').geom('geom1').run('ext1');

model.component('comp1').view('view1').set('renderwireframe', true);
model.component('comp1').view('view1').camera.set('projection', 'orthographic');

model.component('comp1').geom('geom1').run('ext1');
model.component('comp1').geom('geom1').create('mov1', 'Move');
model.component('comp1').geom('geom1').feature('mov1').set('displz', '-3*az_total/2');
model.component('comp1').geom('geom1').feature('mov1').selection('input').set({'ext1'});
model.component('comp1').geom('geom1').run('mov1');
model.component('comp1').geom('geom1').create('arr1', 'Array');
model.component('comp1').geom('geom1').feature('arr1').selection('input').set({'mov1'});
model.component('comp1').geom('geom1').feature('arr1').set('fullsize', [1 1 2]);
model.component('comp1').geom('geom1').feature('arr1').set('displ', {'0' '0' 'az_total'});
model.component('comp1').geom('geom1').runPre('fin');
model.component('comp1').geom('geom1').run;

model.component('comp1').view('view1').set('transparency', true);

model.component('comp1').physics('ewfd').create('pc1', 'PeriodicCondition', 2);
model.component('comp1').physics('ewfd').feature('pc1').selection.set([1 4 7 10 13 16]);
model.component('comp1').physics('ewfd').feature('pc1').set('PeriodicType', 'Floquet');
model.component('comp1').physics('ewfd').feature('pc1').set('kFloquet', [30 40 50]);
model.component('comp1').physics('ewfd').feature('pc1').selection.set([1 4 7 10 13 16 329 330 331 332 333 334]);
model.component('comp1').physics('ewfd').create('pc2', 'PeriodicCondition', 2);
model.component('comp1').physics('ewfd').feature('pc2').selection.set([2 5 8 11 14 17 20 21 22 23 24 25]);
model.component('comp1').physics('ewfd').feature('pc2').set('PeriodicType', 'Floquet');
model.component('comp1').physics('ewfd').feature('pc2').set('kFloquet', [30 40 50]);
model.component('comp1').physics('ewfd').selection.set([1 2 3 4 5 6 7 9 10 12 13 15 16 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 45 46 48 49 51 52 54 55 56 57 58 59 60]);
model.component('comp1').physics('ewfd').feature('wee1').set('DisplacementFieldModel', 'RelativePermittivity');
model.component('comp1').physics('ewfd').feature('wee1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('ewfd').feature('wee1').set('epsilonr', [2.2 0 0 0 2.2 0 0 0 2.2]);
model.component('comp1').physics('ewfd').feature('wee1').set('mur_mat', 'userdef');
model.component('comp1').physics('ewfd').feature('wee1').set('sigma_mat', 'userdef');
model.component('comp1').physics('ewfd').create('pec2', 'PerfectElectricConductor', 2);
model.component('comp1').physics('ewfd').feature('pec2').selection.set([69 78 146 155 177 186 267 276]);

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').selection.set([2 3]);
model.component('comp1').mesh('mesh1').run('edg1');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hauto', 4);
model.component('comp1').mesh('mesh1').run('edg1');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hauto', 2);
model.component('comp1').mesh('mesh1').run('edg1');
model.component('comp1').mesh('mesh1').create('cpe1', 'CopyEdge');
model.component('comp1').mesh('mesh1').feature('cpe1').selection('source').geom(1);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('destination').geom(1);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('source').set([2 3]);
model.component('comp1').mesh('mesh1').feature('cpe1').selection('destination').set([22 507]);
model.component('comp1').mesh('mesh1').run('cpe1');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.set([3 28 47 66 115 140 171 202 221 264]);
model.component('comp1').mesh('mesh1').run('ftri1');
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').create('swe1', 'Sweep');
model.component('comp1').mesh('mesh1').feature('swe1').selection('sourceface').set([3 28 47 66 115 140 171 202 221 264]);
model.component('comp1').mesh('mesh1').feature('swe1').selection('targetface').set([19 44 63 82 131 156 187 218 237 280]);
model.component('comp1').mesh('mesh1').run;
model.component('comp1').mesh('mesh1').feature('swe1').create('dis1', 'Distribution');
model.component('comp1').mesh('mesh1').feature('swe1').feature('dis1').set('numelem', 1);
model.component('comp1').mesh('mesh1').run;

model.study('std1').feature('eig').set('eigunit', 'GHz');
model.study('std1').feature('eig').set('shift', '1');
model.study('std1').feature('eig').set('eigwhich', 'lr');

model.sol.create('sol1');
model.sol('sol1').study('std1');

model.study('std1').feature('eig').set('notlistsolnum', 1);
model.study('std1').feature('eig').set('notsolnum', 'auto');
model.study('std1').feature('eig').set('listsolnum', 1);
model.study('std1').feature('eig').set('solnum', 'auto');

model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'eig');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'eig');
model.sol('sol1').create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').set('control', 'eig');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('e1').create('d1', 'Direct');
model.sol('sol1').feature('e1').feature('d1').set('linsolver', 'mumps');
model.sol('sol1').feature('e1').feature('d1').label([native2unicode(hex2dec({'5e' 'fa'}), 'unicode')  native2unicode(hex2dec({'8b' 'ae'}), 'unicode')  native2unicode(hex2dec({'76' '84'}), 'unicode')  native2unicode(hex2dec({'76' 'f4'}), 'unicode')  native2unicode(hex2dec({'63' 'a5'}), 'unicode')  native2unicode(hex2dec({'6c' '42'}), 'unicode')  native2unicode(hex2dec({'89' 'e3'}), 'unicode')  native2unicode(hex2dec({'56' '68'}), 'unicode') ' (ewfd)']);
model.sol('sol1').attach('std1');

model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').label([native2unicode(hex2dec({'75' '35'}), 'unicode')  native2unicode(hex2dec({'57' '3a'}), 'unicode') ' (ewfd)']);
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').set('data', 'dset1');
model.result('pg1').feature.create('mslc1', 'Multislice');
model.result('pg1').feature('mslc1').set('colortable', 'RainbowLight');
model.result('pg1').feature('mslc1').set('smooth', 'internal');
model.result('pg1').feature('mslc1').set('data', 'parent');
model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical('gev1').label([native2unicode(hex2dec({'72' '79'}), 'unicode')  native2unicode(hex2dec({'5f' '81'}), 'unicode')  native2unicode(hex2dec({'98' '91'}), 'unicode')  native2unicode(hex2dec({'73' '87'}), 'unicode') ' (ewfd)']);
model.result.numerical('gev1').set('data', 'dset1');
model.result.numerical('gev1').set('expr', {'ewfd.freq' 'ewfd.Qfactor'});
model.result.numerical('gev1').set('unit', {'GHz' '1'});
model.result.table.create('tbl1', 'Table');
model.result.numerical('gev1').set('table', 'tbl1');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result.numerical('gev1').setResult;

model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 7);
model.component('comp1').mesh('mesh1').run;

model.sol('sol1').study('std1');

model.study('std1').feature('eig').set('notlistsolnum', 1);
model.study('std1').feature('eig').set('notsolnum', 'auto');
model.study('std1').feature('eig').set('listsolnum', 1);
model.study('std1').feature('eig').set('solnum', 'auto');

model.sol('sol1').feature.remove('e1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'eig');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'eig');
model.sol('sol1').create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').set('control', 'eig');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('e1').create('d1', 'Direct');
model.sol('sol1').feature('e1').feature('d1').set('linsolver', 'mumps');
model.sol('sol1').feature('e1').feature('d1').label([native2unicode(hex2dec({'5e' 'fa'}), 'unicode')  native2unicode(hex2dec({'8b' 'ae'}), 'unicode')  native2unicode(hex2dec({'76' '84'}), 'unicode')  native2unicode(hex2dec({'76' 'f4'}), 'unicode')  native2unicode(hex2dec({'63' 'a5'}), 'unicode')  native2unicode(hex2dec({'6c' '42'}), 'unicode')  native2unicode(hex2dec({'89' 'e3'}), 'unicode')  native2unicode(hex2dec({'56' '68'}), 'unicode') ' (ewfd)']);
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;

out = model;
