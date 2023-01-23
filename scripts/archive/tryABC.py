from abc import ABC, abstractmethod


class Base:
    @staticmethod
    @abstractmethod
    def foo(num):
        pass

    @classmethod
    def bar(cls, nums):
        return [cls.foo(n) for n in nums]

    def fun(self, num):
        return 2*self.foo(num)



class Poop(Base):
    @staticmethod
    def foo(num):
        return num**3

    


nums = list(range(15))
print(Poop.bar(nums))
print(Poop().fun(nums[5]))