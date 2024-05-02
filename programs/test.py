import os
import unittest
import tempfile
import numpy as np
import copy
from termcolor import colored
from show_aj_fct_stellar_params_v2 import probability_composer,query_proba_from_ProductsOdds, query_params_from_ProductsRot, DoClassesExternalvsRot, WriteDataClass

'''
    Test of the probability_composer
'''
class TestProcessProducts(unittest.TestCase):
    def test_probability_composer_basic(self):
        print(colored(" Testing probability_composer()....", "green"), flush=True)
        # Define input
        ProductsOdds = {
                        "header": "# Test header Odds Ratio",
                        "label": ["model1", "model2", "model3", "model4"], 
                        "StarID": ["Star1", "Star2"],
                        "Probabilities": {"key1": np.array([[10, 20, 29, 41], [5, 5, 60, 30]], dtype=float)}, 
                        "Nrows": 2, 
                        "Ncols":3
                        }
        ProductsRot = {"header": "# Test header Rotation",
                       "label": ["model1", "model2", "model3", "model4"],
                        "unit":["U1", "U2", "U3", "U4"],
                        "StarID": ["Star1", "Star2"],
                        "data": {"key0": np.array([[-1, -2, -3, -4],[0,10,20,30]],dtype=float)}
                        }
        Groups = [["model1", "model2"], ["model3"]]
        print(colored("\t Testing with DoRenormalise = False.....", "green"), end="", flush=True)
        RecomposedProductsOdds, RecomposedProductsRot = probability_composer(ProductsOdds, ProductsRot, Groups, DoRenormalise=False)
        # Check results keeping the normalisation over model 1,2,3 and 4
        expected_RecProductsOdds = copy.deepcopy(ProductsOdds)
        expectations= np.array([[30, 30, 29, 41], [10, 10, 60, 30]])
        np.testing.assert_allclose(RecomposedProductsOdds["Probabilities"]["key1"], expectations, atol=1e-8, rtol=1e-5)
        self.assertEqual(RecomposedProductsOdds["header"], expected_RecProductsOdds["header"])
        self.assertEqual(RecomposedProductsOdds["label"], expected_RecProductsOdds["label"])
        print(colored("Passed", "green"), flush=True)

        print(colored("\t Testing with DoRenormalise = True.....", "green"), end="", flush=True)
        RecomposedProductsOdds, RecomposedProductsRot = probability_composer(ProductsOdds, ProductsRot, Groups,  DoRenormalise=True)
        # Check results renormalised on model 1,2 and 3. Model 4 Probability must be 0 after renormalisation in this test
        expected_RecProductsOdds = copy.deepcopy(ProductsOdds)
        norm_star1=30+29
        norm_star2=10+60
        expectations= 100.*np.array([[30/norm_star1, 30/norm_star1, 29/norm_star1, 0.], [10/norm_star2, 10/norm_star2, 60/norm_star2, 0.]])
        np.testing.assert_allclose(RecomposedProductsOdds["Probabilities"]["key1"], expectations, atol=1e-8, rtol=1e-5)
        self.assertEqual(RecomposedProductsOdds["header"], expected_RecProductsOdds["header"])
        self.assertEqual(RecomposedProductsOdds["label"], expected_RecProductsOdds["label"])
        print(colored("Passed", "green"), flush=True)

    def test_probability_composer_second(self):
        # Define the test data
        ProductsOdds = {
            "label": ["1001", "1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"],
            "Probabilities": {
                "key1": np.array([[1.0, 9, 15, 20, 25, 12, 18], [0.5, 4.5, 10, 15, 20, 23, 27]]),
                "key2": np.array([[1, 3, 6, 18, 22, 24, 26], [0.5, 1, 2.5, 5, 7, 16, 68]])
            },
            "Nrows": 2
        }
        ProductsRot = {
            "label": ["1001", "1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"],
            "data": {
                "key1": np.array([[0.1, 1, 2, 3, 4, 5, 6], [2.1, 7, 8, 9, 10, 11, 12]]),
                "key2": np.array([[0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06], [0.01, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12]])
            },
            "Nrows": 2
        }

        expected_RecomposedProductsRot = copy.deepcopy(ProductsRot)

        # Test case 1: Grouping Gate vs Triangle models with DoRenormalise = True and EnforceNorm = False  
        print(colored("\t Second Testing with probability_composer with DoRenormalise = False.....", "green"), end="", flush=True)
            
        Groups = [
            ["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2"],
            ["1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]
        ]
        DoRenormalise = False
        EnforceNorm = False

        expected_RecomposedProductsOdds = {
            "label": ["1001", "1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"],
            "Probabilities": {
                "key1": np.array([[1., 44., 44., 44., 55., 55., 55.], [0.5, 29.5, 29.5, 29.5, 70., 70., 70.]]),
                "key2": np.array([[1., 27., 27., 27., 72., 72., 72.], [0.5, 8.5, 8.5, 8.5, 91., 91., 91.]])
            },
            "Nrows": 2
        }

        RecomposedProductsOdds, RecomposedProductsRot = probability_composer(ProductsOdds, ProductsRot, Groups, DoRenormalise, EnforceNorm)
        
        err_msg = colored("Failed.\n\t Test case 1 failed for RecomposedProductsOdds and DoRenormalise = True: tolerance exceeded in the Probabilities outputs","red")
        np.testing.assert_allclose(RecomposedProductsOdds["Probabilities"]["key1"], expected_RecomposedProductsOdds["Probabilities"]["key1"], atol=1e-8, rtol=1e-5, err_msg=err_msg)
        np.testing.assert_allclose(RecomposedProductsOdds["Probabilities"]["key2"], expected_RecomposedProductsOdds["Probabilities"]["key2"], atol=1e-8, rtol=1e-5, err_msg=err_msg)
        err_msg = colored("Failed.\n\t Test case 1 failed for RecomposedProductsRot and DoRenormalise = True: output data arrays are different","red")
        np.testing.assert_array_equal(RecomposedProductsRot["data"]["key1"], expected_RecomposedProductsRot["data"]["key1"], err_msg=err_msg)
        np.testing.assert_array_equal(RecomposedProductsRot["data"]["key2"], expected_RecomposedProductsRot["data"]["key2"], err_msg=err_msg)
        print(colored("\t Passed", "green"), flush=True)
        
        # --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        # Test case 2: Grouping Gate vs Triangle models with DoRenormalise = True and EnforceNorm = True    
        print(colored("\t Second Testing with probability_composer with DoRenormalise = True.....", "green"), end="", flush=True)
        Groups = [
            ["1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2"],
            ["1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"]
        ]
        DoRenormalise = True
        EnforceNorm = False

        expected_RecomposedProductsOdds = {
            "label": ["1001", "1201Gatedecompose_-1", "1201Gatedecompose_1", "1201Gatedecompose_2", "1201Triangledecompose_-1", "1201Triangledecompose_1", "1201Triangledecompose_2"],
            "Probabilities": {
                "key1": 100.*np.array([[0., 44./99, 44./99, 44./99, 55./99, 55./99, 55./99], [0., 29.5/99.5, 29.5/99.5, 29.5/99.5, 70./99.5, 70./99.5, 70./99.5]]),
                "key2": 100.*np.array([[0., 27./99, 27./99, 27./99, 72./99, 72./99, 72./99], [0., 8.5/99.5, 8.5/99.5, 8.5/99.5, 91./99.5, 91./99.5, 91./99.5]])
            },
            "Nrows": 2
        }

        RecomposedProductsOdds, RecomposedProductsRot = probability_composer(ProductsOdds, ProductsRot, Groups, DoRenormalise, EnforceNorm)

        err_msg = colored("\n\t Test case 2 failed for RecomposedProductsOdds and DoRenormalise = True: tolerance exceeded in the Probabilities outputs","red")
        np.testing.assert_allclose(RecomposedProductsOdds["Probabilities"]["key1"], expected_RecomposedProductsOdds["Probabilities"]["key1"], atol=1e-8, rtol=1e-5, err_msg=err_msg)
        np.testing.assert_allclose(RecomposedProductsOdds["Probabilities"]["key2"], expected_RecomposedProductsOdds["Probabilities"]["key2"], atol=1e-8, rtol=1e-5, err_msg=err_msg)
        err_msg = colored("\n\t Test case 2 failed for RecomposedProductsRot and DoRenormalise = True: output data arrays are different","red")
        np.testing.assert_array_equal(RecomposedProductsRot["data"]["key1"], expected_RecomposedProductsRot["data"]["key1"], err_msg=err_msg)
        np.testing.assert_array_equal(RecomposedProductsRot["data"]["key2"], expected_RecomposedProductsRot["data"]["key2"], err_msg=err_msg)
        print(colored("\t Passed", "green"), flush=True)

    def setUpExternal(self):
        '''
            Defines the test data for External parameters such as Teff
        '''
        self.ExternalData={"starID": ["1435467","3427720","3456181", "3656476", "5184732","6933899", "9139163"], 
                           "label": ["Teff_SDSS","Tot_eTeff","random_eTeff","Teff_IRFM"], 
                           "data": [[6433,86,81,6371], 
                                    [6023,51,34,6303],
                                    [6732,91,65,6554],
                                    [6840,56,32,5786],
                                    [5841,290,287,6137],
                                    [5842,56,36,6014],
                                    [6405,44,39,6623]], 
                           "unit": None, 
                           "Ncols":4., 
                           "Nrows":7}

    def setUpProductsOdds(self):
        '''
            Defines the test data for the Probabilities. Taken from a real file with Tcoef=1
        '''
        self.ProductsOdds = {
            "starID": ["kplr003427720_kasoc-psd_slc_v1", "kplr003656476_kasoc-psd_slc_v1", "kplr009139163_kasoc-wpsd_slc_v1_21000_SK"],
            "label": ["1101", "1111"],
            "Probabilities": {
                "m1s": np.asarray([[45.3 , 48.6], [43.2, 51.5],[34.9, 54.9]], dtype=float),
                "median": np.asarray([[48.2, 51.8], [45.6, 54.4], [38.9, 61.1]],dtype=float),
                "p1s": np.asarray([[51.1, 55.1], [48.1, 57.2],[42.9 , 67.4]], dtype=float)
            }
        }

    def setUpProductsRot(self):
        '''
            Defines the test data for the Rotation parameters. Taken from a real file with Tcoef=1
        '''
        self.ProductsRot = {
            "label": ["a1", "a2", "a2CF"],
            "unit" : ["nHz","nHz","nHz"],
            "starID": ["kplr003427720_kasoc-psd_slc_v1", 
                       "kplr003656476_kasoc-psd_slc_v1", 
                       "kplr003427720_kasoc-psd_slc_v1",
                       "kplr003656476_kasoc-psd_slc_v1",
                       "kplr009139163_kasoc-wpsd_slc_v1_21000_SK",
                       "kplr009139163_kasoc-wpsd_slc_v1_21000_SK"],
            "model": ["1101", "1101", "1111", "1111","1101", "1111"],
            "data": {
                "16": np.asarray([[364.4030, -125.5636, 0], 
                                  [233.1412, -41.0746, 0.], 
                                  [373.8645, -135.7124, 0.], 
                                  [232.9321,-43.2318 ,0.],
                                  [861.7620, -219.2715, 0.0000],
                                  [1025.7865, -220.5088, 0.0000]], dtype=float),
                "50": np.asarray([[533.1651 , -1.7385, 0.], 
                                  [280.8131, 26.0219, 0.], 
                                  [598.4386, -3.4820, 0.], 
                                  [284.3897, 31.9089, 0.],
                                  [1398.9382, -154.7469, 0.0000],
                                  [1823.9419,-154.8882,0.0000]], dtype=float),
                "84": np.asarray([[1199.6619 , 491.7540 , 0.], 
                                  [324.5630  , 77.6375 , 0.], 
                                  [1205.4038 , 498.9627,  0.], 
                                  [342.9008  , 82.2236 , 0.],
                                  [2243.4792 , -29.3571 , 0.0000],
                                  [2913.5499 , -33.5306 , 0.0000]], dtype=float),
            }
        }


    def test_query_proba_from_ProductsOdds_A(self):
        '''
            Querry with ModelCode and a single Confidence provided
        '''
        print(colored("Testing query_params_from_ProductsOdds() with ModelCode and Confidence provided.....", "green"), end="", flush=True)
        self.setUpProductsOdds()
        ID = "kplr003427720_kasoc-psd_slc_v1"
        ModelCode = "1101"
        Confidence = "median"
        result = query_proba_from_ProductsOdds(self.ProductsOdds, ID, ModelCode=ModelCode, Confidence=Confidence)

        expected_result = 48.2
        self.assertEqual(result, expected_result)
        print(colored("Passed", "green"), flush=True)

    def test_query_proba_from_ProductsOdds_B(self):
        '''
            Querry with No ModelCode provided and a single Confidence provided
        '''
        print(colored("Testing query_params_from_ProductsOdds() with No ModelCode and Confidence provided.....", "green"), end="", flush=True)
        self.setUpProductsOdds()
        ID = "kplr003427720_kasoc-psd_slc_v1"
        ModelCode = "1101"
        Confidence = "median"
        result = query_proba_from_ProductsOdds(self.ProductsOdds, ID, Confidence=Confidence)
        self.assertIsInstance(result, np.ndarray, "Result is not a numpy array as expected when no ModelCode is provided")
        expected_result = [48.2, 51.8]
        for i in range(len(expected_result)):
            self.assertEqual(result[i], expected_result[i])
        print(colored("Passed", "green"), flush=True)

    def test_query_proba_from_ProductsOdds_C(self):
        '''
            Querry with No ModelCode provided and No Confidence provided and with an existing ID
        '''
        print(colored("Testing query_params_from_ProductsOdds() with No ModelCode and No Confidence provided.....", "green"), end="", flush=True)
        self.setUpProductsOdds()
        ID = "kplr003427720_kasoc-psd_slc_v1"
        result = query_proba_from_ProductsOdds(self.ProductsOdds, ID)
        self.assertIsInstance(result, dict, "Result is not a dictionary as expected when no optional argument is provided")
        expected_result = {
            "median":[48.2, 51.8],
            "m1s":[45.3,48.6],
            "p1s":[51.1 , 55.1 ],
            "label":["1101", "1111"]
        }
        for key in result.keys():
            for i in range(len(expected_result["label"])):
                self.assertEqual(result[key][i], expected_result[key][i])
        print(colored("Passed", "green"), flush=True)

    def test_query_params_from_ProductsRot_A(self):
        '''
            Querry with all parameters set to existing values, including optional ones
        '''
        print(colored("Testing query_params_from_ProductsRot() with ParamName='a2'.....", "green"), end="", flush=True)

        self.setUpProductsRot()
        ID="kplr003656476_kasoc-psd_slc_v1"
        ModelCode="1101"
        confidence="50"
        ParamName="a2"
        # Test with ResolveID = True
        ResolveID=True
        result = query_params_from_ProductsRot(self.ProductsRot, ID, ModelCode, confidence, ParamName=ParamName, ResolveID=ResolveID)
        expected_result= 26.0219
        self.assertEqual(result, expected_result) 
        # Another test just to be sure
        ModelCode="1111"
        confidence="84"
        result = query_params_from_ProductsRot(self.ProductsRot, ID, ModelCode, confidence, ParamName=ParamName, ResolveID=ResolveID)
        expected_result= 82.2236
        self.assertEqual(result, expected_result)  
        print(colored("Passed", "green"), flush=True)

    def test_query_params_from_ProductsRot_B(self):
        '''
            Querry with all parameters set to existing values, with ParamName=None
        '''
        print(colored("Testing query_params_from_ProductsRot() with ParamName=None.....", "green"), end="", flush=True)

        self.setUpProductsRot()
        ID="kplr003656476_kasoc-psd_slc_v1"
        ModelCode="1101"
        confidence="50"
        ParamName="a2"
        # Test with ResolveID = True
        ResolveID=True
        result = query_params_from_ProductsRot(self.ProductsRot, ID, ModelCode, confidence, ParamName=None, ResolveID=ResolveID)
        self.assertIsInstance(result, np.ndarray, colored("Result is not a numpy array as expected when no ModelCode is provided","red"))
        expected_result= np.asarray([280.8131, 26.0219, 0.], dtype=float)
        for i in range(len(result)):
            self.assertEqual(result[i], expected_result[i]) 
        print(colored("Passed", "green"), flush=True)

    def test_DoClassesExternalvsRot(self):
        print(colored("Testing DoClassesExternalvsRot().....", "green"), end="\n", flush=True)
        self.setUpProductsOdds()
        self.setUpProductsRot()
        self.setUpExternal()
        xlabel_plot="Teff"
        ylabel_plot="a2 (nHz)"
        ModelCode="1111"
        ProbaThresholds=[40, 55]
        # Test 
        ExternalXlabel="Teff_SDSS"
        ExternalXerrlabel="Tot_eTeff"
        RotYlabel="a2"
         # For reference only
        #expected_x_All=[6023,6840,6405]
        #expected_errx_All=[51,56, 44]
        #           KIC                                    P(1101)                   P(1111)
        # kplr003427720_kasoc-psd_slc_v1                     48.2                      51.8                      
        # kplr003656476_kasoc-psd_slc_v1                     45.6                      54.4                      
        # kplr009139163_kasoc-wpsd_slc_v1_21000_SK           38.9                      61.1     
        # a2_16=[-135.7124, -43.2318, -220.5088]
        # a2_50=[-3.4820, 31.9089, -154.8882]
        # a2_84=[498.9627, 82.2236, -33.5306 ]
        # Class for P<ProbaThresholds[0]
        expected_starID_ClassA=[]
        expected_P_ClassA=np.asarray([], dtype=float)
        expected_x_ClassA=np.asarray([] , dtype=float)
        expected_errx_ClassA=np.asarray([[],[]] , dtype=float)
        expected_y_ClassA=np.asarray([] , dtype=float)
        expected_erry_ClassA=np.asarray([[],[]] , dtype=float)
        # Class for ProbaThresholds[0]<P<ProbaThresholds[1]
        expected_starID_ClassB=["003427720", "003656476"]
        expected_P_ClassB=np.asarray([51.8, 54.4], dtype=float)
        expected_x_ClassB=np.asarray([6023, 6840], dtype=float)
        expected_errx_ClassB=np.asarray([[51, 56],[51, 56]], dtype=float)
        expected_y_ClassB=np.asarray([-3.4820 , 31.9089], dtype=float)
        expected_erry_ClassB=np.asarray( [[-3.4820+135.7124, 31.9089+43.2318],
                                          [498.9627+3.4820,  82.2236-31.9089]], dtype=float)
        # Class for P>ProbaThresholds[1]
        expected_starID_ClassC=["009139163"]
        expected_P_ClassC=np.asarray([61.1], dtype=float)
        expected_x_ClassC=np.asarray([6405] , dtype=float)
        expected_errx_ClassC=np.asarray([[44],[44]], dtype=float)
        expected_y_ClassC=np.asarray([-154.8882], dtype=float)
        expected_erry_ClassC=np.asarray( [[-154.8882 + 220.5088],
                                          [-33.5306+154.8882]], dtype=float)
        result = DoClassesExternalvsRot(ModelCode, self.ProductsOdds, self.ProductsRot, self.ExternalData, 
                           RotYlabel, ExternalXlabel, ExternalXerrlabel, xlabel_plot, ylabel_plot,
                           ProbaThresholds=ProbaThresholds)

        self.assertIsInstance(result, dict, "Result is not a dictionary as expected")
        self.assertEqual(result['xlabel'], xlabel_plot)
        self.assertEqual(result['ylabel'], ylabel_plot)
        self.assertTrue('A' in result, "Class 'A' is not found in the result")
        self.assertTrue('B' in result, "Class 'B' is not found in the result")
        self.assertTrue('C' in result, "Class 'C' is not found in the result")       
        #
        self.assertEqual(result['A']['starID'], expected_starID_ClassA)
        self.assertEqual(result['B']['starID'], expected_starID_ClassB)
        self.assertEqual(result['C']['starID'], expected_starID_ClassC)
        #
        self.assertTrue(np.array_equal(result["A"]["Probability"], expected_P_ClassA))
        self.assertTrue(np.array_equal(result["B"]["Probability"], expected_P_ClassB))
        self.assertTrue(np.array_equal(result["C"]["Probability"], expected_P_ClassC))
        #
        self.assertTrue(np.array_equal(result['A']['x'], expected_x_ClassA), "X values do not match for Class 'A'")
        self.assertTrue(np.array_equal(result['A']['err_x'], expected_errx_ClassA), "X error values do not match for Class 'A'")
        self.assertTrue(np.array_equal(result['A']['y'], expected_y_ClassA), "Y values do not match for Class 'A'")
        self.assertTrue(np.array_equal(result['A']['err_y'], expected_erry_ClassA), "Y error values do not match for Class 'A'")
        #
        self.assertTrue(np.array_equal(result['B']['x'], expected_x_ClassB), "X values do not match for Class 'B'")
        self.assertTrue(np.array_equal(result['B']['err_x'], expected_errx_ClassB), "X error values do not match for Class 'B'")
        self.assertTrue(np.array_equal(result['B']['y'], expected_y_ClassB), "Y values do not match for Class 'B'")
        self.assertTrue(np.array_equal(result['B']['err_y'], expected_erry_ClassB), "Y error values do not match for Class 'A'")
        #
        self.assertTrue(np.array_equal(result['C']['x'], expected_x_ClassC), "X values do not match for Class 'C'")
        self.assertTrue(np.array_equal(result['C']['err_x'], expected_errx_ClassC), "X error values do not match for Class 'C'")
        self.assertTrue(np.array_equal(result['C']['y'], expected_y_ClassC), "Y values do not match for Class 'C'")
        self.assertTrue(np.array_equal(result['C']['err_y'], expected_erry_ClassC), "Y error values do not match for Class 'A'")
        print(colored("     Passed", "green"), flush=True)

    def test_WriteDataClass(self):
        print(colored("Testing test_WriteDataClass().....", "green"), flush=True, end="")
        
        temp = tempfile.NamedTemporaryFile(delete=False)
        temp.close()
        data = {
            "key1": {
                "starID": ["star1", "star2"],
                "x": [1.1, 2.2],
                "y": [3.3, 4.4],
                "err_x": [[0.1, 0.2],[0.15,0.25]],
                "err_y": [[0.3, 0.4],[0.35, 0.45]],
                "Probability": [0.9, 0.8],
                "color": 'red',
                "marker": 'o',
                "fillstyle": 'full',
                "label": 'This_is_a_test_label_for_the_unit_test'
            },
            "xlabel": 'X-Label',
            "ylabel": 'Y-Label',
            "hline": True
        }

        WriteDataClass(temp.name, data, header="#This is a test header\n")
        with open(temp.name, 'r') as f:
            content = f.readlines()
        self.assertIn("#This is a test header\n", content[0])
        self.assertIn("x= X-Label\n", content[1])
        self.assertIn("#y= Y-Label\n", content[2])
        labels=content[3][1:].split()
        vals1=content[4].split()
        vals2=content[5].split()
        expected_labels=["key", "StarID", "x", "y", "err_x_inf", "err_x_sup", "err_y_inf", "err_y_sup", "Pr", "color", "marker", "fillstyle", "label"]
        expected_vals1=["key1","star1","1.1","3.3","0.1", "0.15","0.3", "0.35", "0.9","red","o","full","This_is_a_test_label_for_the_unit_test"]
        expected_vals2=["key1", "star2", "2.2", "4.4", "0.2", "0.25", "0.4", "0.45", "0.8", "red", "o", "full", "This_is_a_test_label_for_the_unit_test"]

        for i in range(len(expected_labels)):
            self.assertIn(expected_labels[i], labels[i])
            self.assertIn(expected_vals1[i], vals1[i])
            self.assertIn(expected_vals2[i], vals2[i])
        
        # delete the temporary file
        os.unlink(temp.name)
        print(colored("Passed", "green"), flush=True)

if __name__ == '__main__':
    unittest.main()

