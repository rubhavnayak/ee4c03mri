% FUNCTION TO LOAD THE APPROPRIATE SLICE
% Arguments: Slice number
function [badCh1, badCh2, badCh3, goodCh1, goodCh2, goodCh3] = loadSlice(slice)

    switch slice
        case 1
            badCh1 = load('Slice1/BadData/slice1_channel1.mat', 'slice1_channel1_badData').slice1_channel1_badData;
            badCh2 = load('Slice1/BadData/slice1_channel2.mat', 'slice1_channel2_badData').slice1_channel2_badData;
            badCh3 = load('Slice1/BadData/slice1_channel3.mat', 'slice1_channel3_badData').slice1_channel3_badData;

            goodCh1 = load('Slice1/GoodData/slice1_channel1.mat', 'slice1_channel1_goodData').slice1_channel1_goodData;
            goodCh2 = load('Slice1/GoodData/slice1_channel2.mat', 'slice1_channel2_goodData').slice1_channel2_goodData;
            goodCh3 = load('Slice1/GoodData/slice1_channel3.mat', 'slice1_channel3_goodData').slice1_channel3_goodData;

        case 2
            badCh1 = load('Slice2/BadData/slice2_channel1.mat', 'slice2_channel1_badData').slice2_channel1_badData;
            badCh2 = load('Slice2/BadData/slice2_channel2.mat', 'slice2_channel2_badData').slice2_channel2_badData;
            badCh3 = load('Slice2/BadData/slice2_channel3.mat', 'slice2_channel3_badData').slice2_channel3_badData;

            goodCh1 = load('Slice2/GoodData/slice2_channel1.mat', 'slice2_channel1_goodData').slice2_channel1_goodData;
            goodCh2 = load('Slice2/GoodData/slice2_channel2.mat', 'slice2_channel2_goodData').slice2_channel2_goodData;
            goodCh3 = load('Slice2/GoodData/slice2_channel3.mat', 'slice2_channel3_goodData').slice2_channel3_goodData;

        case 3
            badCh1 = load('Slice3/BadData/slice3_channel1.mat', 'slice3_channel1_badData').slice3_channel1_badData;
            badCh2 = load('Slice3/BadData/slice3_channel2.mat', 'slice3_channel2_badData').slice3_channel2_badData;
            badCh3 = load('Slice3/BadData/slice3_channel3.mat', 'slice3_channel3_badData').slice3_channel3_badData;

            goodCh1 = load('Slice3/GoodData/slice3_channel1.mat', 'slice3_channel1_goodData').slice3_channel1_goodData;
            goodCh2 = load('Slice3/GoodData/slice3_channel2.mat', 'slice3_channel2_goodData').slice3_channel2_goodData;
            goodCh3 = load('Slice3/GoodData/slice3_channel3.mat', 'slice3_channel3_goodData').slice3_channel3_goodData;

        case 4
            badCh1 = load('Slice4/BadData/slice4_channel1.mat', 'slice4_channel1_badData').slice4_channel1_badData;
            badCh2 = load('Slice4/BadData/slice4_channel2.mat', 'slice4_channel2_badData').slice4_channel2_badData;
            badCh3 = load('Slice4/BadData/slice4_channel3.mat', 'slice4_channel3_badData').slice4_channel3_badData;

            goodCh1 = load('Slice4/GoodData/slice4_channel1.mat', 'slice4_channel1_goodData').slice4_channel1_goodData;
            goodCh2 = load('Slice4/GoodData/slice4_channel2.mat', 'slice4_channel2_goodData').slice4_channel2_goodData;
            goodCh3 = load('Slice4/GoodData/slice4_channel3.mat', 'slice4_channel3_goodData').slice4_channel3_goodData;

        case 5
            badCh1 = load('Slice5/BadData/slice5_channel1.mat', 'slice5_channel1_badData').slice5_channel1_badData;
            badCh2 = load('Slice5/BadData/slice5_channel2.mat', 'slice5_channel2_badData').slice5_channel2_badData;
            badCh3 = load('Slice5/BadData/slice5_channel3.mat', 'slice5_channel3_badData').slice5_channel3_badData;

            goodCh1 = load('Slice5/GoodData/slice5_channel1.mat', 'slice5_channel1_goodData').slice5_channel1_goodData;
            goodCh2 = load('Slice5/GoodData/slice5_channel2.mat', 'slice5_channel2_goodData').slice5_channel2_goodData;
            goodCh3 = load('Slice5/GoodData/slice5_channel3.mat', 'slice5_channel3_goodData').slice5_channel3_goodData;

        case 6
            badCh1 = load('Slice6/BadData/slice6_channel1.mat', 'slice6_channel1_badData').slice6_channel1_badData;
            badCh2 = load('Slice6/BadData/slice6_channel2.mat', 'slice6_channel2_badData').slice6_channel2_badData;
            badCh3 = load('Slice6/BadData/slice6_channel3.mat', 'slice6_channel3_badData').slice6_channel3_badData;

            goodCh1 = load('Slice6/GoodData/slice6_channel1.mat', 'slice6_channel1_goodData').slice6_channel1_goodData;
            goodCh2 = load('Slice6/GoodData/slice6_channel2.mat', 'slice6_channel2_goodData').slice6_channel2_goodData;
            goodCh3 = load('Slice6/GoodData/slice6_channel3.mat', 'slice6_channel3_goodData').slice6_channel3_goodData;

        case 7
            badCh1 = load('Slice7/BadData/slice7_channel1.mat', 'slice7_channel1_badData').slice7_channel1_badData;
            badCh2 = load('Slice7/BadData/slice7_channel2.mat', 'slice7_channel2_badData').slice7_channel2_badData;
            badCh3 = load('Slice7/BadData/slice7_channel3.mat', 'slice7_channel3_badData').slice7_channel3_badData;

            goodCh1 = load('Slice7/GoodData/slice7_channel1.mat', 'slice7_channel1_goodData').slice7_channel1_goodData;
            goodCh2 = load('Slice7/GoodData/slice7_channel2.mat', 'slice7_channel2_goodData').slice7_channel2_goodData;
            goodCh3 = load('Slice7/GoodData/slice7_channel3.mat', 'slice7_channel3_goodData').slice7_channel3_goodData;

        case 8
            badCh1 = load('Slice8/BadData/slice8_channel1.mat', 'slice8_channel1_badData').slice8_channel1_badData;
            badCh2 = load('Slice8/BadData/slice8_channel2.mat', 'slice8_channel2_badData').slice8_channel2_badData;
            badCh3 = load('Slice8/BadData/slice8_channel3.mat', 'slice8_channel3_badData').slice8_channel3_badData;

            goodCh1 = load('Slice8/GoodData/slice8_channel1.mat', 'slice8_channel1_goodData').slice8_channel1_goodData;
            goodCh2 = load('Slice8/GoodData/slice8_channel2.mat', 'slice8_channel2_goodData').slice8_channel2_goodData;
            goodCh3 = load('Slice8/GoodData/slice8_channel3.mat', 'slice8_channel3_goodData').slice8_channel3_goodData;

        otherwise
            error('INVALID SLICE');
    end
end