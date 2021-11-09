--- test.py	(original)
+++ test.py	(refactored)
@@ -24,7 +24,7 @@
 
 
 if __name__ == "__main__":
-    list_of_numbers = range(0, 5)
+    list_of_numbers = list(range(0, 5))
     list_of_objects = [MyClass(i) for i in list_of_numbers]
 
     pool = mp.Pool(NUM_CORE)
@@ -33,5 +33,5 @@
     pool.close()
     pool.join()
 
-    print list_of_numbers
-    print list_of_results
+    print(list_of_numbers)
+    print(list_of_results)
