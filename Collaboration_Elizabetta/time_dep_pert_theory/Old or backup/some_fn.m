function tmp = some_fn(in1,in2,in3,in4)
%SOME_FN
%    TMP = SOME_FN(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 5.11.
%    15-Oct-2014 14:13:47

Fge1_1 = in1(1,:);
Fge1_2 = in1(5,:);
Fge1_3 = in1(9,:);
Fge1_4 = in1(13,:);
Fge1_5 = in1(17,:);
Fge1_6 = in1(21,:);
Fge1_7 = in1(25,:);
Fge1_8 = in1(29,:);
Fge2_1 = in1(2,:);
Fge2_2 = in1(6,:);
Fge2_3 = in1(10,:);
Fge2_4 = in1(14,:);
Fge2_5 = in1(18,:);
Fge2_6 = in1(22,:);
Fge2_7 = in1(26,:);
Fge2_8 = in1(30,:);
Fge3_1 = in1(3,:);
Fge3_2 = in1(7,:);
Fge3_3 = in1(11,:);
Fge3_4 = in1(15,:);
Fge3_5 = in1(19,:);
Fge3_6 = in1(23,:);
Fge3_7 = in1(27,:);
Fge3_8 = in1(31,:);
Fge4_1 = in1(4,:);
Fge4_2 = in1(8,:);
Fge4_3 = in1(12,:);
Fge4_4 = in1(16,:);
Fge4_5 = in1(20,:);
Fge4_6 = in1(24,:);
Fge4_7 = in1(28,:);
Fge4_8 = in1(32,:);
Gge1_1 = in2(1,:);
Gge1_2 = in2(5,:);
Gge1_3 = in2(9,:);
Gge1_4 = in2(13,:);
Gge1_5 = in2(17,:);
Gge1_6 = in2(21,:);
Gge1_7 = in2(25,:);
Gge1_8 = in2(29,:);
Gge2_1 = in2(2,:);
Gge2_2 = in2(6,:);
Gge2_3 = in2(10,:);
Gge2_4 = in2(14,:);
Gge2_5 = in2(18,:);
Gge2_6 = in2(22,:);
Gge2_7 = in2(26,:);
Gge2_8 = in2(30,:);
Gge3_1 = in2(3,:);
Gge3_2 = in2(7,:);
Gge3_3 = in2(11,:);
Gge3_4 = in2(15,:);
Gge3_5 = in2(19,:);
Gge3_6 = in2(23,:);
Gge3_7 = in2(27,:);
Gge3_8 = in2(31,:);
Gge4_1 = in2(4,:);
Gge4_2 = in2(8,:);
Gge4_3 = in2(12,:);
Gge4_4 = in2(16,:);
Gge4_5 = in2(20,:);
Gge4_6 = in2(24,:);
Gge4_7 = in2(28,:);
Gge4_8 = in2(32,:);
Hge1_1 = in3(1,:);
Hge1_2 = in3(5,:);
Hge1_3 = in3(9,:);
Hge1_4 = in3(13,:);
Hge1_5 = in3(17,:);
Hge1_6 = in3(21,:);
Hge1_7 = in3(25,:);
Hge1_8 = in3(29,:);
Hge2_1 = in3(2,:);
Hge2_2 = in3(6,:);
Hge2_3 = in3(10,:);
Hge2_4 = in3(14,:);
Hge2_5 = in3(18,:);
Hge2_6 = in3(22,:);
Hge2_7 = in3(26,:);
Hge2_8 = in3(30,:);
Hge3_1 = in3(3,:);
Hge3_2 = in3(7,:);
Hge3_3 = in3(11,:);
Hge3_4 = in3(15,:);
Hge3_5 = in3(19,:);
Hge3_6 = in3(23,:);
Hge3_7 = in3(27,:);
Hge3_8 = in3(31,:);
Hge4_1 = in3(4,:);
Hge4_2 = in3(8,:);
Hge4_3 = in3(12,:);
Hge4_4 = in3(16,:);
Hge4_5 = in3(20,:);
Hge4_6 = in3(24,:);
Hge4_7 = in3(28,:);
Hge4_8 = in3(32,:);
Ige1_1 = in4(1,:);
Ige1_2 = in4(5,:);
Ige1_3 = in4(9,:);
Ige1_4 = in4(13,:);
Ige1_5 = in4(17,:);
Ige1_6 = in4(21,:);
Ige1_7 = in4(25,:);
Ige1_8 = in4(29,:);
Ige2_1 = in4(2,:);
Ige2_2 = in4(6,:);
Ige2_3 = in4(10,:);
Ige2_4 = in4(14,:);
Ige2_5 = in4(18,:);
Ige2_6 = in4(22,:);
Ige2_7 = in4(26,:);
Ige2_8 = in4(30,:);
Ige3_1 = in4(3,:);
Ige3_2 = in4(7,:);
Ige3_3 = in4(11,:);
Ige3_4 = in4(15,:);
Ige3_5 = in4(19,:);
Ige3_6 = in4(23,:);
Ige3_7 = in4(27,:);
Ige3_8 = in4(31,:);
Ige4_1 = in4(4,:);
Ige4_2 = in4(8,:);
Ige4_3 = in4(12,:);
Ige4_4 = in4(16,:);
Ige4_5 = in4(20,:);
Ige4_6 = in4(24,:);
Ige4_7 = in4(28,:);
Ige4_8 = in4(32,:);
t2 = conj(Hge1_1);
t3 = Fge1_1.*t2;
t4 = conj(Hge1_2);
t5 = Fge1_2.*t4;
t6 = conj(Hge1_3);
t7 = Fge1_3.*t6;
t8 = conj(Hge1_4);
t9 = Fge1_4.*t8;
t10 = conj(Hge1_5);
t11 = Fge1_5.*t10;
t12 = conj(Hge1_6);
t13 = Fge1_6.*t12;
t14 = conj(Hge1_7);
t15 = Fge1_7.*t14;
t16 = conj(Hge1_8);
t17 = Fge1_8.*t16;
t18 = t3+t5+t7+t9+t11+t13+t15+t17;
t19 = conj(Hge2_1);
t20 = Fge1_1.*t19;
t21 = conj(Hge2_2);
t22 = Fge1_2.*t21;
t23 = conj(Hge2_3);
t24 = Fge1_3.*t23;
t25 = conj(Hge2_4);
t26 = Fge1_4.*t25;
t27 = conj(Hge2_5);
t28 = Fge1_5.*t27;
t29 = conj(Hge2_6);
t30 = Fge1_6.*t29;
t31 = conj(Hge2_7);
t32 = Fge1_7.*t31;
t33 = conj(Hge2_8);
t34 = Fge1_8.*t33;
t35 = t20+t22+t24+t26+t28+t30+t32+t34;
t36 = conj(Hge3_1);
t37 = Fge1_1.*t36;
t38 = conj(Hge3_2);
t39 = Fge1_2.*t38;
t40 = conj(Hge3_3);
t41 = Fge1_3.*t40;
t42 = conj(Hge3_4);
t43 = Fge1_4.*t42;
t44 = conj(Hge3_5);
t45 = Fge1_5.*t44;
t46 = conj(Hge3_6);
t47 = Fge1_6.*t46;
t48 = conj(Hge3_7);
t49 = Fge1_7.*t48;
t50 = conj(Hge3_8);
t51 = Fge1_8.*t50;
t52 = t37+t39+t41+t43+t45+t47+t49+t51;
t53 = conj(Hge4_1);
t54 = Fge1_1.*t53;
t55 = conj(Hge4_2);
t56 = Fge1_2.*t55;
t57 = conj(Hge4_3);
t58 = Fge1_3.*t57;
t59 = conj(Hge4_4);
t60 = Fge1_4.*t59;
t61 = conj(Hge4_5);
t62 = Fge1_5.*t61;
t63 = conj(Hge4_6);
t64 = Fge1_6.*t63;
t65 = conj(Hge4_7);
t66 = Fge1_7.*t65;
t67 = conj(Hge4_8);
t68 = Fge1_8.*t67;
t69 = t54+t56+t58+t60+t62+t64+t66+t68;
t70 = Fge2_1.*t2;
t71 = Fge2_2.*t4;
t72 = Fge2_3.*t6;
t73 = Fge2_4.*t8;
t74 = Fge2_5.*t10;
t75 = Fge2_6.*t12;
t76 = Fge2_7.*t14;
t77 = Fge2_8.*t16;
t78 = t70+t71+t72+t73+t74+t75+t76+t77;
t79 = Fge2_1.*t19;
t80 = Fge2_2.*t21;
t81 = Fge2_3.*t23;
t82 = Fge2_4.*t25;
t83 = Fge2_5.*t27;
t84 = Fge2_6.*t29;
t85 = Fge2_7.*t31;
t86 = Fge2_8.*t33;
t87 = t79+t80+t81+t82+t83+t84+t85+t86;
t88 = Fge2_1.*t36;
t89 = Fge2_2.*t38;
t90 = Fge2_3.*t40;
t91 = Fge2_4.*t42;
t92 = Fge2_5.*t44;
t93 = Fge2_6.*t46;
t94 = Fge2_7.*t48;
t95 = Fge2_8.*t50;
t96 = t88+t89+t90+t91+t92+t93+t94+t95;
t97 = Fge2_1.*t53;
t98 = Fge2_2.*t55;
t99 = Fge2_3.*t57;
t100 = Fge2_4.*t59;
t101 = Fge2_5.*t61;
t102 = Fge2_6.*t63;
t103 = Fge2_7.*t65;
t104 = Fge2_8.*t67;
t105 = t97+t98+t99+t100+t101+t102+t103+t104;
t106 = Fge3_1.*t2;
t107 = Fge3_2.*t4;
t108 = Fge3_3.*t6;
t109 = Fge3_4.*t8;
t110 = Fge3_5.*t10;
t111 = Fge3_6.*t12;
t112 = Fge3_7.*t14;
t113 = Fge3_8.*t16;
t114 = t106+t107+t108+t109+t110+t111+t112+t113;
t115 = Fge3_1.*t19;
t116 = Fge3_2.*t21;
t117 = Fge3_3.*t23;
t118 = Fge3_4.*t25;
t119 = Fge3_5.*t27;
t120 = Fge3_6.*t29;
t121 = Fge3_7.*t31;
t122 = Fge3_8.*t33;
t123 = t115+t116+t117+t118+t119+t120+t121+t122;
t124 = Fge3_1.*t36;
t125 = Fge3_2.*t38;
t126 = Fge3_3.*t40;
t127 = Fge3_4.*t42;
t128 = Fge3_5.*t44;
t129 = Fge3_6.*t46;
t130 = Fge3_7.*t48;
t131 = Fge3_8.*t50;
t132 = t124+t125+t126+t127+t128+t129+t130+t131;
t133 = Fge3_1.*t53;
t134 = Fge3_2.*t55;
t135 = Fge3_3.*t57;
t136 = Fge3_4.*t59;
t137 = Fge3_5.*t61;
t138 = Fge3_6.*t63;
t139 = Fge3_7.*t65;
t140 = Fge3_8.*t67;
t141 = t133+t134+t135+t136+t137+t138+t139+t140;
t142 = Fge4_1.*t2;
t143 = Fge4_2.*t4;
t144 = Fge4_3.*t6;
t145 = Fge4_4.*t8;
t146 = Fge4_5.*t10;
t147 = Fge4_6.*t12;
t148 = Fge4_7.*t14;
t149 = Fge4_8.*t16;
t150 = t142+t143+t144+t145+t146+t147+t148+t149;
t151 = Fge4_1.*t19;
t152 = Fge4_2.*t21;
t153 = Fge4_3.*t23;
t154 = Fge4_4.*t25;
t155 = Fge4_5.*t27;
t156 = Fge4_6.*t29;
t157 = Fge4_7.*t31;
t158 = Fge4_8.*t33;
t159 = t151+t152+t153+t154+t155+t156+t157+t158;
t160 = Fge4_1.*t36;
t161 = Fge4_2.*t38;
t162 = Fge4_3.*t40;
t163 = Fge4_4.*t42;
t164 = Fge4_5.*t44;
t165 = Fge4_6.*t46;
t166 = Fge4_7.*t48;
t167 = Fge4_8.*t50;
t168 = t160+t161+t162+t163+t164+t165+t166+t167;
t169 = Fge4_1.*t53;
t170 = Fge4_2.*t55;
t171 = Fge4_3.*t57;
t172 = Fge4_4.*t59;
t173 = Fge4_5.*t61;
t174 = Fge4_6.*t63;
t175 = Fge4_7.*t65;
t176 = Fge4_8.*t67;
t177 = t169+t170+t171+t172+t173+t174+t175+t176;
t178 = Gge1_1.*t2;
t179 = Gge1_2.*t4;
t180 = Gge1_3.*t6;
t181 = Gge1_4.*t8;
t182 = Gge1_5.*t10;
t183 = Gge1_6.*t12;
t184 = Gge1_7.*t14;
t185 = Gge1_8.*t16;
t186 = t178+t179+t180+t181+t182+t183+t184+t185;
t187 = Gge1_1.*t19;
t188 = Gge1_2.*t21;
t189 = Gge1_3.*t23;
t190 = Gge1_4.*t25;
t191 = Gge1_5.*t27;
t192 = Gge1_6.*t29;
t193 = Gge1_7.*t31;
t194 = Gge1_8.*t33;
t195 = t187+t188+t189+t190+t191+t192+t193+t194;
t196 = Gge1_1.*t36;
t197 = Gge1_2.*t38;
t198 = Gge1_3.*t40;
t199 = Gge1_4.*t42;
t200 = Gge1_5.*t44;
t201 = Gge1_6.*t46;
t202 = Gge1_7.*t48;
t203 = Gge1_8.*t50;
t204 = t196+t197+t198+t199+t200+t201+t202+t203;
t205 = Gge1_1.*t53;
t206 = Gge1_2.*t55;
t207 = Gge1_3.*t57;
t208 = Gge1_4.*t59;
t209 = Gge1_5.*t61;
t210 = Gge1_6.*t63;
t211 = Gge1_7.*t65;
t212 = Gge1_8.*t67;
t213 = t205+t206+t207+t208+t209+t210+t211+t212;
t214 = Gge2_1.*t2;
t215 = Gge2_2.*t4;
t216 = Gge2_3.*t6;
t217 = Gge2_4.*t8;
t218 = Gge2_5.*t10;
t219 = Gge2_6.*t12;
t220 = Gge2_7.*t14;
t221 = Gge2_8.*t16;
t222 = t214+t215+t216+t217+t218+t219+t220+t221;
t223 = Gge2_1.*t19;
t224 = Gge2_2.*t21;
t225 = Gge2_3.*t23;
t226 = Gge2_4.*t25;
t227 = Gge2_5.*t27;
t228 = Gge2_6.*t29;
t229 = Gge2_7.*t31;
t230 = Gge2_8.*t33;
t231 = t223+t224+t225+t226+t227+t228+t229+t230;
t232 = Gge2_1.*t36;
t233 = Gge2_2.*t38;
t234 = Gge2_3.*t40;
t235 = Gge2_4.*t42;
t236 = Gge2_5.*t44;
t237 = Gge2_6.*t46;
t238 = Gge2_7.*t48;
t239 = Gge2_8.*t50;
t240 = t232+t233+t234+t235+t236+t237+t238+t239;
t241 = Gge2_1.*t53;
t242 = Gge2_2.*t55;
t243 = Gge2_3.*t57;
t244 = Gge2_4.*t59;
t245 = Gge2_5.*t61;
t246 = Gge2_6.*t63;
t247 = Gge2_7.*t65;
t248 = Gge2_8.*t67;
t249 = t241+t242+t243+t244+t245+t246+t247+t248;
t250 = Gge3_1.*t2;
t251 = Gge3_2.*t4;
t252 = Gge3_3.*t6;
t253 = Gge3_4.*t8;
t254 = Gge3_5.*t10;
t255 = Gge3_6.*t12;
t256 = Gge3_7.*t14;
t257 = Gge3_8.*t16;
t258 = t250+t251+t252+t253+t254+t255+t256+t257;
t259 = Gge3_1.*t19;
t260 = Gge3_2.*t21;
t261 = Gge3_3.*t23;
t262 = Gge3_4.*t25;
t263 = Gge3_5.*t27;
t264 = Gge3_6.*t29;
t265 = Gge3_7.*t31;
t266 = Gge3_8.*t33;
t267 = t259+t260+t261+t262+t263+t264+t265+t266;
t268 = Gge3_1.*t36;
t269 = Gge3_2.*t38;
t270 = Gge3_3.*t40;
t271 = Gge3_4.*t42;
t272 = Gge3_5.*t44;
t273 = Gge3_6.*t46;
t274 = Gge3_7.*t48;
t275 = Gge3_8.*t50;
t276 = t268+t269+t270+t271+t272+t273+t274+t275;
t277 = Gge3_1.*t53;
t278 = Gge3_2.*t55;
t279 = Gge3_3.*t57;
t280 = Gge3_4.*t59;
t281 = Gge3_5.*t61;
t282 = Gge3_6.*t63;
t283 = Gge3_7.*t65;
t284 = Gge3_8.*t67;
t285 = t277+t278+t279+t280+t281+t282+t283+t284;
t286 = Gge4_1.*t2;
t287 = Gge4_2.*t4;
t288 = Gge4_3.*t6;
t289 = Gge4_4.*t8;
t290 = Gge4_5.*t10;
t291 = Gge4_6.*t12;
t292 = Gge4_7.*t14;
t293 = Gge4_8.*t16;
t294 = t286+t287+t288+t289+t290+t291+t292+t293;
t295 = Gge4_1.*t19;
t296 = Gge4_2.*t21;
t297 = Gge4_3.*t23;
t298 = Gge4_4.*t25;
t299 = Gge4_5.*t27;
t300 = Gge4_6.*t29;
t301 = Gge4_7.*t31;
t302 = Gge4_8.*t33;
t303 = t295+t296+t297+t298+t299+t300+t301+t302;
t304 = Gge4_1.*t36;
t305 = Gge4_2.*t38;
t306 = Gge4_3.*t40;
t307 = Gge4_4.*t42;
t308 = Gge4_5.*t44;
t309 = Gge4_6.*t46;
t310 = Gge4_7.*t48;
t311 = Gge4_8.*t50;
t312 = t304+t305+t306+t307+t308+t309+t310+t311;
t313 = Gge4_1.*t53;
t314 = Gge4_2.*t55;
t315 = Gge4_3.*t57;
t316 = Gge4_4.*t59;
t317 = Gge4_5.*t61;
t318 = Gge4_6.*t63;
t319 = Gge4_7.*t65;
t320 = Gge4_8.*t67;
t321 = t313+t314+t315+t316+t317+t318+t319+t320;
t322 = conj(Fge1_1);
t323 = conj(Fge1_2);
t324 = Ige1_1.*t2;
t325 = Ige1_2.*t4;
t326 = Ige1_3.*t6;
t327 = Ige1_4.*t8;
t328 = Ige1_5.*t10;
t329 = Ige1_6.*t12;
t330 = Ige1_7.*t14;
t331 = Ige1_8.*t16;
t332 = t324+t325+t326+t327+t328+t329+t330+t331;
t333 = Ige1_1.*t19;
t334 = Ige1_2.*t21;
t335 = Ige1_3.*t23;
t336 = Ige1_4.*t25;
t337 = Ige1_5.*t27;
t338 = Ige1_6.*t29;
t339 = Ige1_7.*t31;
t340 = Ige1_8.*t33;
t341 = t333+t334+t335+t336+t337+t338+t339+t340;
t342 = Ige1_1.*t36;
t343 = Ige1_2.*t38;
t344 = Ige1_3.*t40;
t345 = Ige1_4.*t42;
t346 = Ige1_5.*t44;
t347 = Ige1_6.*t46;
t348 = Ige1_7.*t48;
t349 = Ige1_8.*t50;
t350 = t342+t343+t344+t345+t346+t347+t348+t349;
t351 = Ige1_1.*t53;
t352 = Ige1_2.*t55;
t353 = Ige1_3.*t57;
t354 = Ige1_4.*t59;
t355 = Ige1_5.*t61;
t356 = Ige1_6.*t63;
t357 = Ige1_7.*t65;
t358 = Ige1_8.*t67;
t359 = t351+t352+t353+t354+t355+t356+t357+t358;
t360 = conj(Fge1_3);
t361 = conj(Fge1_4);
t362 = conj(Fge1_5);
t363 = conj(Fge1_6);
t364 = conj(Fge1_7);
t365 = conj(Fge1_8);
t366 = conj(Fge2_1);
t367 = conj(Fge2_2);
t368 = Ige2_1.*t2;
t369 = Ige2_2.*t4;
t370 = Ige2_3.*t6;
t371 = Ige2_4.*t8;
t372 = Ige2_5.*t10;
t373 = Ige2_6.*t12;
t374 = Ige2_7.*t14;
t375 = Ige2_8.*t16;
t376 = t368+t369+t370+t371+t372+t373+t374+t375;
t377 = Ige2_1.*t19;
t378 = Ige2_2.*t21;
t379 = Ige2_3.*t23;
t380 = Ige2_4.*t25;
t381 = Ige2_5.*t27;
t382 = Ige2_6.*t29;
t383 = Ige2_7.*t31;
t384 = Ige2_8.*t33;
t385 = t377+t378+t379+t380+t381+t382+t383+t384;
t386 = Ige2_1.*t36;
t387 = Ige2_2.*t38;
t388 = Ige2_3.*t40;
t389 = Ige2_4.*t42;
t390 = Ige2_5.*t44;
t391 = Ige2_6.*t46;
t392 = Ige2_7.*t48;
t393 = Ige2_8.*t50;
t394 = t386+t387+t388+t389+t390+t391+t392+t393;
t395 = Ige2_1.*t53;
t396 = Ige2_2.*t55;
t397 = Ige2_3.*t57;
t398 = Ige2_4.*t59;
t399 = Ige2_5.*t61;
t400 = Ige2_6.*t63;
t401 = Ige2_7.*t65;
t402 = Ige2_8.*t67;
t403 = t395+t396+t397+t398+t399+t400+t401+t402;
t404 = conj(Fge2_3);
t405 = conj(Fge2_4);
t406 = conj(Fge2_5);
t407 = conj(Fge2_6);
t408 = conj(Fge2_7);
t409 = conj(Fge2_8);
t410 = conj(Fge3_1);
t411 = conj(Fge3_2);
t412 = Ige3_1.*t2;
t413 = Ige3_2.*t4;
t414 = Ige3_3.*t6;
t415 = Ige3_4.*t8;
t416 = Ige3_5.*t10;
t417 = Ige3_6.*t12;
t418 = Ige3_7.*t14;
t419 = Ige3_8.*t16;
t420 = t412+t413+t414+t415+t416+t417+t418+t419;
t421 = Ige3_1.*t19;
t422 = Ige3_2.*t21;
t423 = Ige3_3.*t23;
t424 = Ige3_4.*t25;
t425 = Ige3_5.*t27;
t426 = Ige3_6.*t29;
t427 = Ige3_7.*t31;
t428 = Ige3_8.*t33;
t429 = t421+t422+t423+t424+t425+t426+t427+t428;
t430 = Ige3_1.*t36;
t431 = Ige3_2.*t38;
t432 = Ige3_3.*t40;
t433 = Ige3_4.*t42;
t434 = Ige3_5.*t44;
t435 = Ige3_6.*t46;
t436 = Ige3_7.*t48;
t437 = Ige3_8.*t50;
t438 = t430+t431+t432+t433+t434+t435+t436+t437;
t439 = Ige3_1.*t53;
t440 = Ige3_2.*t55;
t441 = Ige3_3.*t57;
t442 = Ige3_4.*t59;
t443 = Ige3_5.*t61;
t444 = Ige3_6.*t63;
t445 = Ige3_7.*t65;
t446 = Ige3_8.*t67;
t447 = t439+t440+t441+t442+t443+t444+t445+t446;
t448 = conj(Fge3_3);
t449 = conj(Fge3_4);
t450 = conj(Fge3_5);
t451 = conj(Fge3_6);
t452 = conj(Fge3_7);
t453 = conj(Fge3_8);
t454 = conj(Fge4_1);
t455 = conj(Fge4_2);
t456 = Ige4_1.*t2;
t457 = Ige4_2.*t4;
t458 = Ige4_3.*t6;
t459 = Ige4_4.*t8;
t460 = Ige4_5.*t10;
t461 = Ige4_6.*t12;
t462 = Ige4_7.*t14;
t463 = Ige4_8.*t16;
t464 = t456+t457+t458+t459+t460+t461+t462+t463;
t465 = Ige4_1.*t19;
t466 = Ige4_2.*t21;
t467 = Ige4_3.*t23;
t468 = Ige4_4.*t25;
t469 = Ige4_5.*t27;
t470 = Ige4_6.*t29;
t471 = Ige4_7.*t31;
t472 = Ige4_8.*t33;
t473 = t465+t466+t467+t468+t469+t470+t471+t472;
t474 = Ige4_1.*t36;
t475 = Ige4_2.*t38;
t476 = Ige4_3.*t40;
t477 = Ige4_4.*t42;
t478 = Ige4_5.*t44;
t479 = Ige4_6.*t46;
t480 = Ige4_7.*t48;
t481 = Ige4_8.*t50;
t482 = t474+t475+t476+t477+t478+t479+t480+t481;
t483 = Ige4_1.*t53;
t484 = Ige4_2.*t55;
t485 = Ige4_3.*t57;
t486 = Ige4_4.*t59;
t487 = Ige4_5.*t61;
t488 = Ige4_6.*t63;
t489 = Ige4_7.*t65;
t490 = Ige4_8.*t67;
t491 = t483+t484+t485+t486+t487+t488+t489+t490;
t492 = conj(Fge4_3);
t493 = conj(Fge4_4);
t494 = conj(Fge4_5);
t495 = conj(Fge4_6);
t496 = conj(Fge4_7);
t497 = conj(Fge4_8);
t498 = conj(Gge1_1);
t499 = conj(Gge1_2);
t500 = conj(Gge1_3);
t501 = conj(Gge1_4);
t502 = conj(Gge1_5);
t503 = conj(Gge1_6);
t504 = conj(Gge1_7);
t505 = conj(Gge1_8);
t506 = conj(Gge2_1);
t507 = conj(Gge2_2);
t508 = conj(Gge2_3);
t509 = conj(Gge2_4);
t510 = conj(Gge2_5);
t511 = conj(Gge2_6);
t512 = conj(Gge2_7);
t513 = conj(Gge2_8);
t514 = conj(Gge3_1);
t515 = conj(Gge3_2);
t516 = conj(Gge3_3);
t517 = conj(Gge3_4);
t518 = conj(Gge3_5);
t519 = conj(Gge3_6);
t520 = conj(Gge3_7);
t521 = conj(Gge3_8);
t522 = conj(Gge4_1);
t523 = conj(Gge4_2);
t524 = conj(Gge4_3);
t525 = conj(Gge4_4);
t526 = conj(Gge4_5);
t527 = conj(Gge4_6);
t528 = conj(Gge4_7);
t529 = conj(Gge4_8);
t530 = Fge1_1.*t498;
t531 = Fge1_2.*t499;
t532 = Fge1_3.*t500;
t533 = Fge1_4.*t501;
t534 = Fge1_5.*t502;
t535 = Fge1_6.*t503;
t536 = Fge1_7.*t504;
t537 = Fge1_8.*t505;
t538 = t530+t531+t532+t533+t534+t535+t536+t537;
t539 = Fge1_1.*t506;
t540 = Fge1_2.*t507;
t541 = Fge1_3.*t508;
t542 = Fge1_4.*t509;
t543 = Fge1_5.*t510;
t544 = Fge1_6.*t511;
t545 = Fge1_7.*t512;
t546 = Fge1_8.*t513;
t547 = t539+t540+t541+t542+t543+t544+t545+t546;
t548 = Fge1_1.*t514;
t549 = Fge1_2.*t515;
t550 = Fge1_3.*t516;
t551 = Fge1_4.*t517;
t552 = Fge1_5.*t518;
t553 = Fge1_6.*t519;
t554 = Fge1_7.*t520;
t555 = Fge1_8.*t521;
t556 = t548+t549+t550+t551+t552+t553+t554+t555;
t557 = Fge1_1.*t522;
t558 = Fge1_2.*t523;
t559 = Fge1_3.*t524;
t560 = Fge1_4.*t525;
t561 = Fge1_5.*t526;
t562 = Fge1_6.*t527;
t563 = Fge1_7.*t528;
t564 = Fge1_8.*t529;
t565 = t557+t558+t559+t560+t561+t562+t563+t564;
t566 = Fge2_1.*t498;
t567 = Fge2_2.*t499;
t568 = Fge2_3.*t500;
t569 = Fge2_4.*t501;
t570 = Fge2_5.*t502;
t571 = Fge2_6.*t503;
t572 = Fge2_7.*t504;
t573 = Fge2_8.*t505;
t574 = t566+t567+t568+t569+t570+t571+t572+t573;
t575 = Fge2_1.*t506;
t576 = Fge2_2.*t507;
t577 = Fge2_3.*t508;
t578 = Fge2_4.*t509;
t579 = Fge2_5.*t510;
t580 = Fge2_6.*t511;
t581 = Fge2_7.*t512;
t582 = Fge2_8.*t513;
t583 = t575+t576+t577+t578+t579+t580+t581+t582;
t584 = Fge2_1.*t514;
t585 = Fge2_2.*t515;
t586 = Fge2_3.*t516;
t587 = Fge2_4.*t517;
t588 = Fge2_5.*t518;
t589 = Fge2_6.*t519;
t590 = Fge2_7.*t520;
t591 = Fge2_8.*t521;
t592 = t584+t585+t586+t587+t588+t589+t590+t591;
t593 = Fge2_1.*t522;
t594 = Fge2_2.*t523;
t595 = Fge2_3.*t524;
t596 = Fge2_4.*t525;
t597 = Fge2_5.*t526;
t598 = Fge2_6.*t527;
t599 = Fge2_7.*t528;
t600 = Fge2_8.*t529;
t601 = t593+t594+t595+t596+t597+t598+t599+t600;
t602 = Fge3_1.*t498;
t603 = Fge3_2.*t499;
t604 = Fge3_3.*t500;
t605 = Fge3_4.*t501;
t606 = Fge3_5.*t502;
t607 = Fge3_6.*t503;
t608 = Fge3_7.*t504;
t609 = Fge3_8.*t505;
t610 = t602+t603+t604+t605+t606+t607+t608+t609;
t611 = Fge3_1.*t506;
t612 = Fge3_2.*t507;
t613 = Fge3_3.*t508;
t614 = Fge3_4.*t509;
t615 = Fge3_5.*t510;
t616 = Fge3_6.*t511;
t617 = Fge3_7.*t512;
t618 = Fge3_8.*t513;
t619 = t611+t612+t613+t614+t615+t616+t617+t618;
t620 = Fge3_1.*t514;
t621 = Fge3_2.*t515;
t622 = Fge3_3.*t516;
t623 = Fge3_4.*t517;
t624 = Fge3_5.*t518;
t625 = Fge3_6.*t519;
t626 = Fge3_7.*t520;
t627 = Fge3_8.*t521;
t628 = t620+t621+t622+t623+t624+t625+t626+t627;
t629 = Fge3_1.*t522;
t630 = Fge3_2.*t523;
t631 = Fge3_3.*t524;
t632 = Fge3_4.*t525;
t633 = Fge3_5.*t526;
t634 = Fge3_6.*t527;
t635 = Fge3_7.*t528;
t636 = Fge3_8.*t529;
t637 = t629+t630+t631+t632+t633+t634+t635+t636;
t638 = Fge4_1.*t498;
t639 = Fge4_2.*t499;
t640 = Fge4_3.*t500;
t641 = Fge4_4.*t501;
t642 = Fge4_5.*t502;
t643 = Fge4_6.*t503;
t644 = Fge4_7.*t504;
t645 = Fge4_8.*t505;
t646 = t638+t639+t640+t641+t642+t643+t644+t645;
t647 = Fge4_1.*t506;
t648 = Fge4_2.*t507;
t649 = Fge4_3.*t508;
t650 = Fge4_4.*t509;
t651 = Fge4_5.*t510;
t652 = Fge4_6.*t511;
t653 = Fge4_7.*t512;
t654 = Fge4_8.*t513;
t655 = t647+t648+t649+t650+t651+t652+t653+t654;
t656 = Fge4_1.*t514;
t657 = Fge4_2.*t515;
t658 = Fge4_3.*t516;
t659 = Fge4_4.*t517;
t660 = Fge4_5.*t518;
t661 = Fge4_6.*t519;
t662 = Fge4_7.*t520;
t663 = Fge4_8.*t521;
t664 = t656+t657+t658+t659+t660+t661+t662+t663;
t665 = Fge4_1.*t522;
t666 = Fge4_2.*t523;
t667 = Fge4_3.*t524;
t668 = Fge4_4.*t525;
t669 = Fge4_5.*t526;
t670 = Fge4_6.*t527;
t671 = Fge4_7.*t528;
t672 = Fge4_8.*t529;
t673 = t665+t666+t667+t668+t669+t670+t671+t672;
tmp = imag(t322.*(Gge1_1.*t332+Gge2_1.*t341+Gge3_1.*t350+Gge4_1.*t359)).*1.961488763309136+imag(t323.*(Gge1_2.*t332+Gge2_2.*t341+Gge3_2.*t350+Gge4_2.*t359)).*1.961488763309136+imag(t360.*(Gge1_3.*t332+Gge2_3.*t341+Gge3_3.*t350+Gge4_3.*t359)).*1.961488763309136+imag(t361.*(Gge1_4.*t332+Gge2_4.*t341+Gge3_4.*t350+Gge4_4.*t359)).*1.961488763309136+imag(t362.*(Gge1_5.*t332+Gge2_5.*t341+Gge3_5.*t350+Gge4_5.*t359)).*1.961488763309136+imag(t363.*(Gge1_6.*t332+Gge2_6.*t341+Gge3_6.*t350+Gge4_6.*t359)).*1.961488763309136+imag(t364.*(Gge1_7.*t332+Gge2_7.*t341+Gge3_7.*t350+Gge4_7.*t359)).*1.961488763309136+imag(t365.*(Gge1_8.*t332+Gge2_8.*t341+Gge3_8.*t350+Gge4_8.*t359)).*1.961488763309136+imag(t366.*(Gge1_1.*t376+Gge2_1.*t385+Gge3_1.*t394+Gge4_1.*t403)).*1.358033504315331e-2+imag(t367.*(Gge1_2.*t376+Gge2_2.*t385+Gge3_2.*t394+Gge4_2.*t403)).*1.358033504315331e-2+imag(t404.*(Gge1_3.*t376+Gge2_3.*t385+Gge3_3.*t394+Gge4_3.*t403)).*1.358033504315331e-2+imag(t405.*(Gge1_4.*t376+Gge2_4.*t385+Gge3_4.*t394+Gge4_4.*t403)).*1.358033504315331e-2+imag(t406.*(Gge1_5.*t376+Gge2_5.*t385+Gge3_5.*t394+Gge4_5.*t403)).*1.358033504315331e-2+imag(t407.*(Gge1_6.*t376+Gge2_6.*t385+Gge3_6.*t394+Gge4_6.*t403)).*1.358033504315331e-2+imag(t408.*(Gge1_7.*t376+Gge2_7.*t385+Gge3_7.*t394+Gge4_7.*t403)).*1.358033504315331e-2+imag(t409.*(Gge1_8.*t376+Gge2_8.*t385+Gge3_8.*t394+Gge4_8.*t403)).*1.358033504315331e-2+imag(t410.*(Gge1_1.*t420+Gge2_1.*t429+Gge3_1.*t438+Gge4_1.*t447)).*1.358033504315332e-2+imag(t411.*(Gge1_2.*t420+Gge2_2.*t429+Gge3_2.*t438+Gge4_2.*t447)).*1.358033504315332e-2+imag(t448.*(Gge1_3.*t420+Gge2_3.*t429+Gge3_3.*t438+Gge4_3.*t447)).*1.358033504315332e-2+imag(t449.*(Gge1_4.*t420+Gge2_4.*t429+Gge3_4.*t438+Gge4_4.*t447)).*1.358033504315332e-2+imag(t450.*(Gge1_5.*t420+Gge2_5.*t429+Gge3_5.*t438+Gge4_5.*t447)).*1.358033504315332e-2+imag(t451.*(Gge1_6.*t420+Gge2_6.*t429+Gge3_6.*t438+Gge4_6.*t447)).*1.358033504315332e-2+imag(t452.*(Gge1_7.*t420+Gge2_7.*t429+Gge3_7.*t438+Gge4_7.*t447)).*1.358033504315332e-2+imag(t453.*(Gge1_8.*t420+Gge2_8.*t429+Gge3_8.*t438+Gge4_8.*t447)).*1.358033504315332e-2+imag(t454.*(Gge1_1.*t464+Gge2_1.*t473+Gge3_1.*t482+Gge4_1.*t491)).*9.408400100109731e-5+imag(t455.*(Gge1_2.*t464+Gge2_2.*t473+Gge3_2.*t482+Gge4_2.*t491)).*9.408400100109731e-5+imag(t492.*(Gge1_3.*t464+Gge2_3.*t473+Gge3_3.*t482+Gge4_3.*t491)).*9.408400100109731e-5+imag(t493.*(Gge1_4.*t464+Gge2_4.*t473+Gge3_4.*t482+Gge4_4.*t491)).*9.408400100109731e-5+imag(t494.*(Gge1_5.*t464+Gge2_5.*t473+Gge3_5.*t482+Gge4_5.*t491)).*9.408400100109731e-5+imag(t495.*(Gge1_6.*t464+Gge2_6.*t473+Gge3_6.*t482+Gge4_6.*t491)).*9.408400100109731e-5+imag(t496.*(Gge1_7.*t464+Gge2_7.*t473+Gge3_7.*t482+Gge4_7.*t491)).*9.408400100109731e-5+imag(t497.*(Gge1_8.*t464+Gge2_8.*t473+Gge3_8.*t482+Gge4_8.*t491)).*9.408400100109731e-5+imag(t498.*(Ige1_1.*t18+Ige2_1.*t35+Ige3_1.*t52+Ige4_1.*t69)).*1.961488763309136+imag(t499.*(Ige1_2.*t18+Ige2_2.*t35+Ige3_2.*t52+Ige4_2.*t69)).*1.961488763309136+imag(t500.*(Ige1_3.*t18+Ige2_3.*t35+Ige3_3.*t52+Ige4_3.*t69)).*1.961488763309136+imag(t501.*(Ige1_4.*t18+Ige2_4.*t35+Ige3_4.*t52+Ige4_4.*t69)).*1.961488763309136+imag(t502.*(Ige1_5.*t18+Ige2_5.*t35+Ige3_5.*t52+Ige4_5.*t69)).*1.961488763309136+imag(t503.*(Ige1_6.*t18+Ige2_6.*t35+Ige3_6.*t52+Ige4_6.*t69)).*1.961488763309136+imag(t504.*(Ige1_7.*t18+Ige2_7.*t35+Ige3_7.*t52+Ige4_7.*t69)).*1.961488763309136+imag(t505.*(Ige1_8.*t18+Ige2_8.*t35+Ige3_8.*t52+Ige4_8.*t69)).*1.961488763309136+imag(t506.*(Ige1_1.*t78+Ige2_1.*t87+Ige3_1.*t96+Ige4_1.*t105)).*1.358033504315331e-2+imag(t507.*(Ige1_2.*t78+Ige2_2.*t87+Ige3_2.*t96+Ige4_2.*t105)).*1.358033504315331e-2+imag(t508.*(Ige1_3.*t78+Ige2_3.*t87+Ige3_3.*t96+Ige4_3.*t105)).*1.358033504315331e-2+imag(t509.*(Ige1_4.*t78+Ige2_4.*t87+Ige3_4.*t96+Ige4_4.*t105)).*1.358033504315331e-2+imag(t510.*(Ige1_5.*t78+Ige2_5.*t87+Ige3_5.*t96+Ige4_5.*t105)).*1.358033504315331e-2+imag(t511.*(Ige1_6.*t78+Ige2_6.*t87+Ige3_6.*t96+Ige4_6.*t105)).*1.358033504315331e-2+imag(t512.*(Ige1_7.*t78+Ige2_7.*t87+Ige3_7.*t96+Ige4_7.*t105)).*1.358033504315331e-2+imag(t513.*(Ige1_8.*t78+Ige2_8.*t87+Ige3_8.*t96+Ige4_8.*t105)).*1.358033504315331e-2+imag(t514.*(Ige1_1.*t114+Ige2_1.*t123+Ige3_1.*t132+Ige4_1.*t141)).*1.358033504315332e-2+imag(t515.*(Ige1_2.*t114+Ige2_2.*t123+Ige3_2.*t132+Ige4_2.*t141)).*1.358033504315332e-2+imag(t516.*(Ige1_3.*t114+Ige2_3.*t123+Ige3_3.*t132+Ige4_3.*t141)).*1.358033504315332e-2+imag(t517.*(Ige1_4.*t114+Ige2_4.*t123+Ige3_4.*t132+Ige4_4.*t141)).*1.358033504315332e-2+imag(t518.*(Ige1_5.*t114+Ige2_5.*t123+Ige3_5.*t132+Ige4_5.*t141)).*1.358033504315332e-2+imag(t519.*(Ige1_6.*t114+Ige2_6.*t123+Ige3_6.*t132+Ige4_6.*t141)).*1.358033504315332e-2+imag(t520.*(Ige1_7.*t114+Ige2_7.*t123+Ige3_7.*t132+Ige4_7.*t141)).*1.358033504315332e-2+imag(t521.*(Ige1_8.*t114+Ige2_8.*t123+Ige3_8.*t132+Ige4_8.*t141)).*1.358033504315332e-2+imag(t322.*(Ige1_1.*t186+Ige2_1.*t195+Ige3_1.*t204+Ige4_1.*t213)).*1.961488763309136+imag(t323.*(Ige1_2.*t186+Ige2_2.*t195+Ige3_2.*t204+Ige4_2.*t213)).*1.961488763309136+imag(t360.*(Ige1_3.*t186+Ige2_3.*t195+Ige3_3.*t204+Ige4_3.*t213)).*1.961488763309136+imag(t361.*(Ige1_4.*t186+Ige2_4.*t195+Ige3_4.*t204+Ige4_4.*t213)).*1.961488763309136+imag(t362.*(Ige1_5.*t186+Ige2_5.*t195+Ige3_5.*t204+Ige4_5.*t213)).*1.961488763309136+imag(t522.*(Ige1_1.*t150+Ige2_1.*t159+Ige3_1.*t168+Ige4_1.*t177)).*9.408400100109731e-5+imag(t363.*(Ige1_6.*t186+Ige2_6.*t195+Ige3_6.*t204+Ige4_6.*t213)).*1.961488763309136+imag(t523.*(Ige1_2.*t150+Ige2_2.*t159+Ige3_2.*t168+Ige4_2.*t177)).*9.408400100109731e-5+imag(t364.*(Ige1_7.*t186+Ige2_7.*t195+Ige3_7.*t204+Ige4_7.*t213)).*1.961488763309136+imag(t524.*(Ige1_3.*t150+Ige2_3.*t159+Ige3_3.*t168+Ige4_3.*t177)).*9.408400100109731e-5+imag(t365.*(Ige1_8.*t186+Ige2_8.*t195+Ige3_8.*t204+Ige4_8.*t213)).*1.961488763309136+imag(t525.*(Ige1_4.*t150+Ige2_4.*t159+Ige3_4.*t168+Ige4_4.*t177)).*9.408400100109731e-5+imag(t526.*(Ige1_5.*t150+Ige2_5.*t159+Ige3_5.*t168+Ige4_5.*t177)).*9.408400100109731e-5+imag(t527.*(Ige1_6.*t150+Ige2_6.*t159+Ige3_6.*t168+Ige4_6.*t177)).*9.408400100109731e-5+imag(t528.*(Ige1_7.*t150+Ige2_7.*t159+Ige3_7.*t168+Ige4_7.*t177)).*9.408400100109731e-5+imag(t529.*(Ige1_8.*t150+Ige2_8.*t159+Ige3_8.*t168+Ige4_8.*t177)).*9.408400100109731e-5+imag(t366.*(Ige1_1.*t222+Ige2_1.*t231+Ige3_1.*t240+Ige4_1.*t249)).*1.358033504315331e-2+imag(t367.*(Ige1_2.*t222+Ige2_2.*t231+Ige3_2.*t240+Ige4_2.*t249)).*1.358033504315331e-2+imag(t404.*(Ige1_3.*t222+Ige2_3.*t231+Ige3_3.*t240+Ige4_3.*t249)).*1.358033504315331e-2+imag(t405.*(Ige1_4.*t222+Ige2_4.*t231+Ige3_4.*t240+Ige4_4.*t249)).*1.358033504315331e-2+imag(t406.*(Ige1_5.*t222+Ige2_5.*t231+Ige3_5.*t240+Ige4_5.*t249)).*1.358033504315331e-2+imag(t407.*(Ige1_6.*t222+Ige2_6.*t231+Ige3_6.*t240+Ige4_6.*t249)).*1.358033504315331e-2+imag(t408.*(Ige1_7.*t222+Ige2_7.*t231+Ige3_7.*t240+Ige4_7.*t249)).*1.358033504315331e-2+imag(t409.*(Ige1_8.*t222+Ige2_8.*t231+Ige3_8.*t240+Ige4_8.*t249)).*1.358033504315331e-2+imag(t410.*(Ige1_1.*t258+Ige2_1.*t267+Ige3_1.*t276+Ige4_1.*t285)).*1.358033504315332e-2+imag(t411.*(Ige1_2.*t258+Ige2_2.*t267+Ige3_2.*t276+Ige4_2.*t285)).*1.358033504315332e-2+imag(t448.*(Ige1_3.*t258+Ige2_3.*t267+Ige3_3.*t276+Ige4_3.*t285)).*1.358033504315332e-2+imag(t449.*(Ige1_4.*t258+Ige2_4.*t267+Ige3_4.*t276+Ige4_4.*t285)).*1.358033504315332e-2+imag(t450.*(Ige1_5.*t258+Ige2_5.*t267+Ige3_5.*t276+Ige4_5.*t285)).*1.358033504315332e-2+imag(t451.*(Ige1_6.*t258+Ige2_6.*t267+Ige3_6.*t276+Ige4_6.*t285)).*1.358033504315332e-2+imag(t452.*(Ige1_7.*t258+Ige2_7.*t267+Ige3_7.*t276+Ige4_7.*t285)).*1.358033504315332e-2+imag(t453.*(Ige1_8.*t258+Ige2_8.*t267+Ige3_8.*t276+Ige4_8.*t285)).*1.358033504315332e-2+imag(t454.*(Ige1_1.*t294+Ige2_1.*t303+Ige3_1.*t312+Ige4_1.*t321)).*9.408400100109731e-5+imag(t455.*(Ige1_2.*t294+Ige2_2.*t303+Ige3_2.*t312+Ige4_2.*t321)).*9.408400100109731e-5+imag(t492.*(Ige1_3.*t294+Ige2_3.*t303+Ige3_3.*t312+Ige4_3.*t321)).*9.408400100109731e-5+imag(t493.*(Ige1_4.*t294+Ige2_4.*t303+Ige3_4.*t312+Ige4_4.*t321)).*9.408400100109731e-5+imag(t494.*(Ige1_5.*t294+Ige2_5.*t303+Ige3_5.*t312+Ige4_5.*t321)).*9.408400100109731e-5+imag(t495.*(Ige1_6.*t294+Ige2_6.*t303+Ige3_6.*t312+Ige4_6.*t321)).*9.408400100109731e-5+imag(t496.*(Ige1_7.*t294+Ige2_7.*t303+Ige3_7.*t312+Ige4_7.*t321)).*9.408400100109731e-5+imag(t497.*(Ige1_8.*t294+Ige2_8.*t303+Ige3_8.*t312+Ige4_8.*t321)).*9.408400100109731e-5+imag(t2.*(Ige1_1.*t538+Ige2_1.*t547+Ige3_1.*t556+Ige4_1.*t565)).*1.961488763309136+imag(t4.*(Ige1_2.*t538+Ige2_2.*t547+Ige3_2.*t556+Ige4_2.*t565)).*1.961488763309136+imag(t6.*(Ige1_3.*t538+Ige2_3.*t547+Ige3_3.*t556+Ige4_3.*t565)).*1.961488763309136+imag(t8.*(Ige1_4.*t538+Ige2_4.*t547+Ige3_4.*t556+Ige4_4.*t565)).*1.961488763309136+imag(t10.*(Ige1_5.*t538+Ige2_5.*t547+Ige3_5.*t556+Ige4_5.*t565)).*1.961488763309136+imag(t12.*(Ige1_6.*t538+Ige2_6.*t547+Ige3_6.*t556+Ige4_6.*t565)).*1.961488763309136+imag(t14.*(Ige1_7.*t538+Ige2_7.*t547+Ige3_7.*t556+Ige4_7.*t565)).*1.961488763309136+imag(t16.*(Ige1_8.*t538+Ige2_8.*t547+Ige3_8.*t556+Ige4_8.*t565)).*1.961488763309136+imag(t19.*(Ige1_1.*t574+Ige2_1.*t583+Ige3_1.*t592+Ige4_1.*t601)).*1.358033504315331e-2+imag(t21.*(Ige1_2.*t574+Ige2_2.*t583+Ige3_2.*t592+Ige4_2.*t601)).*1.358033504315331e-2+imag(t23.*(Ige1_3.*t574+Ige2_3.*t583+Ige3_3.*t592+Ige4_3.*t601)).*1.358033504315331e-2+imag(t25.*(Ige1_4.*t574+Ige2_4.*t583+Ige3_4.*t592+Ige4_4.*t601)).*1.358033504315331e-2+imag(t27.*(Ige1_5.*t574+Ige2_5.*t583+Ige3_5.*t592+Ige4_5.*t601)).*1.358033504315331e-2+imag(t29.*(Ige1_6.*t574+Ige2_6.*t583+Ige3_6.*t592+Ige4_6.*t601)).*1.358033504315331e-2+imag(t31.*(Ige1_7.*t574+Ige2_7.*t583+Ige3_7.*t592+Ige4_7.*t601)).*1.358033504315331e-2+imag(t33.*(Ige1_8.*t574+Ige2_8.*t583+Ige3_8.*t592+Ige4_8.*t601)).*1.358033504315331e-2+imag(t36.*(Ige1_1.*t610+Ige2_1.*t619+Ige3_1.*t628+Ige4_1.*t637)).*1.358033504315332e-2+imag(t38.*(Ige1_2.*t610+Ige2_2.*t619+Ige3_2.*t628+Ige4_2.*t637)).*1.358033504315332e-2+imag(t40.*(Ige1_3.*t610+Ige2_3.*t619+Ige3_3.*t628+Ige4_3.*t637)).*1.358033504315332e-2+imag(t42.*(Ige1_4.*t610+Ige2_4.*t619+Ige3_4.*t628+Ige4_4.*t637)).*1.358033504315332e-2+imag(t44.*(Ige1_5.*t610+Ige2_5.*t619+Ige3_5.*t628+Ige4_5.*t637)).*1.358033504315332e-2+imag(t46.*(Ige1_6.*t610+Ige2_6.*t619+Ige3_6.*t628+Ige4_6.*t637)).*1.358033504315332e-2+imag(t48.*(Ige1_7.*t610+Ige2_7.*t619+Ige3_7.*t628+Ige4_7.*t637)).*1.358033504315332e-2+imag(t50.*(Ige1_8.*t610+Ige2_8.*t619+Ige3_8.*t628+Ige4_8.*t637)).*1.358033504315332e-2+imag(t53.*(Ige1_1.*t646+Ige2_1.*t655+Ige3_1.*t664+Ige4_1.*t673)).*9.408400100109731e-5+imag(t55.*(Ige1_2.*t646+Ige2_2.*t655+Ige3_2.*t664+Ige4_2.*t673)).*9.408400100109731e-5+imag(t57.*(Ige1_3.*t646+Ige2_3.*t655+Ige3_3.*t664+Ige4_3.*t673)).*9.408400100109731e-5+imag(t59.*(Ige1_4.*t646+Ige2_4.*t655+Ige3_4.*t664+Ige4_4.*t673)).*9.408400100109731e-5+imag(t61.*(Ige1_5.*t646+Ige2_5.*t655+Ige3_5.*t664+Ige4_5.*t673)).*9.408400100109731e-5+imag(t63.*(Ige1_6.*t646+Ige2_6.*t655+Ige3_6.*t664+Ige4_6.*t673)).*9.408400100109731e-5+imag(t65.*(Ige1_7.*t646+Ige2_7.*t655+Ige3_7.*t664+Ige4_7.*t673)).*9.408400100109731e-5+imag(t67.*(Ige1_8.*t646+Ige2_8.*t655+Ige3_8.*t664+Ige4_8.*t673)).*9.408400100109731e-5;