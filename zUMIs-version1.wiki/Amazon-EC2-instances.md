![](https://d0.awsstatic.com/logos/aws/Powered-by-Amazon-Web-Services(!).png)


If you are a user of the Amazon Cloud (Amazon Web Services), we are providing a machine image (AMI) with a complete, functional zUMIs installation.

You can start your EC2 instance with zUMIs in the **US East (N. Virginia)** location (change this in the top right corner).
If you need to use a different location, the AMI may be copied through the AWS console.

In the **AWS Console**, click on **AMIs**.
Here, switch the search bar to **Public images** and find the zUMIs image by pasting 
> ** ami-5c98fe23** (AMI ID) 

or 

> ** zUMIs_25052018** (image name).

Klick the **Launch** button to start the instance. When selecting the instance type, remember STAR needs at least 30 Gb RAM for loading the genome index. 
Follow the assistant to launch the instance and you are good to go!

To use the new instance, connect to it through ssh. You will get instructions from Amazon if you click on the **"Connect"** button in the EC2 Dashboard. It should look like this:
```
ssh -i "your-key.pem" ec2-user@ec2-your-instance.compute-1.amazonaws.com
```


Once logged in, you can find the zUMIs installation in the following path of the EC2 instance:
```
~/programs/zUMIs/
```

The example dataset and its output are located in:
```
~/data/example/
```