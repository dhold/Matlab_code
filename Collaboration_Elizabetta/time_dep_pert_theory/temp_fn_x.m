function out =temp_fn_x(om,om1,om2,om_u,t_sep) 
out =exp(t_sep.*(om1 + om2).*(0.0 - 1.0.*i)).*((0.07938596823200929 - 0.1019696817481795.*i)./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(- 94.89557798119951.*i - 17956.36774133398) + om2.*(- 55.06384876029245.*i - 17332.06223846845) + om1.*om2 + om1.^2 + (724742.6737125896.*i + 10818308.54317013))) + (0.000807589109646587 - 0.1708561387403608.*i)./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(104.2630681233358.*i + 16083.45123273737) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10433310.88222013 - 625274.0027475406.*i))) + (0.000950597352336443.*i - 0.00007802728056820856)./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(18580.67324419952 - 104.2630681233358.*i) + om2.*(17956.36774133398 - 64.43133890242872.*i) + om1.*om2 + om1.^2 + (11207692.7807475 - 755458.0170976262.*i))) + (0.1019696817481795.*i - 0.07938596823200929)./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(85.52808783906324.*i + 17332.06223846845) + om2.*(55.06384876029245.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1677.478253230826 - 528008.087760838.*i))) + (0.000807589109646587 - 0.002221929074551924.*i)./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(17332.06223846845 - 122.9980484076083.*i) + om2.*(17332.06223846845 - 73.79882904456498.*i) + om1.*om2 + om1.^2 - (852723.9320843287.*i + 3630.844778899285))) + (0.1708561387403608.*i + 0.000807589109646587)./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(16083.45123273737 - 104.2630681233358.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (625274.0027475406.*i + 10433310.88222013))) + (- 0.000950597352336443.*i - 0.00007802728056820856)./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(- 104.2630681233358.*i - 18580.67324419952) + om2.*(- 64.43133890242872.*i - 17956.36774133398) + om1.*om2 + om1.^2 + (755458.0170976262.*i + 11207692.7807475))) + (83.98854743614871.*i - 0.07938596823200929)./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(76.16059769692698.*i + 16707.75673560291) + om2.*(45.69635861815619.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1392.104793972758 - 508989.0956633501.*i))) + (- 0.1708561387403608.*i - 0.000807589109646587)./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(16707.75673560291 - 113.630558265472.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (822008.5886992921.*i + 3169.971576515177))) + (0.0785003518417945 - 0.2737764178408768.*i)./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(15459.14572987184 - 94.89557798119951.*i) + om2.*(16083.45123273737 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (606255.0106500527.*i + 10043180.39798112))) + (0.1708561387403608.*i - 0.000807589109646587)./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(113.630558265472.*i + 16707.75673560291) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3169.971576515177 - 822008.5886992921.*i))) + (- 83.98854743614871.*i - 0.07938596823200929)./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(16707.75673560291 - 76.16059769692698.*i) + om2.*(16707.75673560291 - 45.69635861815619.*i) + om1.*om2 + om1.^2 - (508989.0956633501.*i + 1392.104793972758))) + (0.1019696817481795.*i + 0.07938596823200929)./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(17332.06223846845 - 85.52808783906324.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (528008.087760838.*i + 1677.478253230826))) + (0.002221929074551924.*i + 0.000807589109646587)./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(122.9980484076083.*i + 17332.06223846845) + om2.*(73.79882904456498.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3630.844778899285 - 852723.9320843287.*i))) + (0.2737764178408768.*i + 0.0785003518417945)./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(94.89557798119951.*i + 15459.14572987184) + om2.*(55.06384876029245.*i + 16083.45123273737) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10043180.39798112 - 606255.0106500527.*i))) + (- 0.1019696817481795.*i - 0.07938596823200929)./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(17956.36774133398 - 94.89557798119951.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 + (10818308.54317013 - 724742.6737125896.*i))) + (0.1708561387403608.*i - 0.000807589109646587)./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(104.2630681233358.*i + 16083.45123273737) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10433310.88222013 - 625274.0027475406.*i))) + (0.002221929074551924.*i - 0.000807589109646587)./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(17332.06223846845 - 122.9980484076083.*i) + om2.*(17332.06223846845 - 73.79882904456498.*i) + om1.*om2 + om1.^2 - (852723.9320843287.*i + 3630.844778899285))) + (0.07938596823200929 - 0.1019696817481795.*i)./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(85.52808783906324.*i + 17332.06223846845) + om2.*(55.06384876029245.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1677.478253230826 - 528008.087760838.*i))) + (0.1019696817481795.*i - 0.07938596823200929)./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(- 94.89557798119951.*i - 17956.36774133398) + om2.*(- 55.06384876029245.*i - 17332.06223846845) + om1.*om2 + om1.^2 + (724742.6737125896.*i + 10818308.54317013))) + (0.00007802728056820856 - 0.000950597352336443.*i)./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(18580.67324419952 - 104.2630681233358.*i) + om2.*(17956.36774133398 - 64.43133890242872.*i) + om1.*om2 + om1.^2 + (11207692.7807475 - 755458.0170976262.*i))) + (0.07938596823200929 - 83.98854743614871.*i)./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(76.16059769692698.*i + 16707.75673560291) + om2.*(45.69635861815619.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (1392.104793972758 - 508989.0956633501.*i))) + (0.1708561387403608.*i + 0.000807589109646587)./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(16707.75673560291 - 113.630558265472.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (822008.5886992921.*i + 3169.971576515177))) + (0.000807589109646587 - 0.1708561387403608.*i)./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(113.630558265472.*i + 16707.75673560291) + om2.*(64.43133890242872.*i + 16707.75673560291) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3169.971576515177 - 822008.5886992921.*i))) + (- 0.1708561387403608.*i - 0.000807589109646587)./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(16083.45123273737 - 104.2630681233358.*i) + om2.*(16707.75673560291 - 64.43133890242872.*i) + om1.*om2 + om1.^2 - (625274.0027475406.*i + 10433310.88222013))) + (0.000950597352336443.*i + 0.00007802728056820856)./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(- 104.2630681233358.*i - 18580.67324419952) + om2.*(- 64.43133890242872.*i - 17956.36774133398) + om1.*om2 + om1.^2 + (755458.0170976262.*i + 11207692.7807475))) + (83.98854743614871.*i + 0.07938596823200929)./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(16707.75673560291 - 76.16059769692698.*i) + om2.*(16707.75673560291 - 45.69635861815619.*i) + om1.*om2 + om1.^2 - (508989.0956633501.*i + 1392.104793972758))) + (- 0.2737764178408768.*i - 0.0785003518417945)./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(94.89557798119951.*i + 15459.14572987184) + om2.*(55.06384876029245.*i + 16083.45123273737) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (10043180.39798112 - 606255.0106500527.*i))) + (0.2737764178408768.*i - 0.0785003518417945)./((om.*(1.0.*i + 0.0) + (16707.75673560291.*i + 15.2321195393854)).*(om1.*(15459.14572987184 - 94.89557798119951.*i) + om2.*(16083.45123273737 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (606255.0106500527.*i + 10043180.39798112))) + (- 0.1019696817481795.*i - 0.07938596823200929)./((om.*(1.0.*i + 0.0) + (15.2321195393854 - 16707.75673560291.*i)).*(om1.*(17332.06223846845 - 85.52808783906324.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 - (528008.087760838.*i + 1677.478253230826))) + (- 0.002221929074551924.*i - 0.000807589109646587)./((om.*(1.0.*i + 0.0) + (24.59960968152166 - 17332.06223846845.*i)).*(om1.*(122.9980484076083.*i + 17332.06223846845) + om2.*(73.79882904456498.*i + 17332.06223846845) - 1.0.*om1.*om2 - 1.0.*om1.^2 + (3630.844778899285 - 852723.9320843287.*i))) + (0.1019696817481795.*i + 0.07938596823200929)./((om.*(1.0.*i + 0.0) + (17332.06223846845.*i + 24.59960968152166)).*(om1.*(17956.36774133398 - 94.89557798119951.*i) + om2.*(17332.06223846845 - 55.06384876029245.*i) + om1.*om2 + om1.^2 + (10818308.54317013 - 724742.6737125896.*i))));