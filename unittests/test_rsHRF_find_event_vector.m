function tests = test_rsHRF_find_event_vector
% Unit Tests for rsHRF_find_event_vector
tests = functiontests(localfunctions);

function test_rsHRF_find_event_vector_1(testCase)
import matlab.unittest.constraints.*
matrix = randn(200,1);
ts = [20,60,100 150];
matrix(ts)=10;  % plot(matrix)
event = rsHRF_find_event_vector(matrix,1,1,[]);
testCase.verifyThat(event, IsOfClass('double'));
testCase.assertTrue(all(ismember(ts,find(event))));
